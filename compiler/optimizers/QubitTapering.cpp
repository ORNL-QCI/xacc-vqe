#include "QubitTapering.hpp"
#include "PauliOperator.hpp"

#include <algorithm>

namespace xacc {
namespace vqe {
unsigned pivot(Eigen::MatrixXd &B, unsigned r, unsigned c) {
  for (unsigned k = r; k < B.rows(); k++) {
    if (std::fabs(B(k, c)) > 1e-12) {
      return k;
    }
  }
  return B.rows();
}

void reduced_row_echelon_form(const Eigen::MatrixXd &A, Eigen::MatrixXd &B,
                              std::vector<int> &pivotCols) {
  unsigned row = A.rows(), col = A.cols();
  unsigned index = 0, i, j, k;
  B = A;

  for (i = 0; i < col; i++) {
    if (index == row)
      break;

    k = pivot(B, index, i);
    if (k == row)
      continue;
    if (k != index) {
      B.row(k).swap(B.row(index));
    }

    B.row(index) *= 1.0 / B(index, i);

    for (j = 0; j < row; j++) {
      if (j == index)
        continue;

      B.row(j) += (-1.0 * B(j, i)) * B.row(index);
    }

    index++;
  }

  unsigned row2 = 0;
  for (unsigned col = 0; col < B.cols() && row2 < B.rows(); col++) {
    if (std::fabs(B(row2, col)) < 1e-12) 
      continue;

    pivotCols.push_back(col);
    row2++;
  }
}

std::shared_ptr<IR> QubitTapering::transform(std::shared_ptr<IR> ir) {

  PauliOperator H;
  H.fromXACCIR(ir);

  auto n = H.nQubits();

  int nRows = 0;
  Eigen::MatrixXi tableaux;
  for (auto &term : H) {
    Eigen::VectorXi x(n), z(n), row(2 * n);
    auto p = term.second.toBinaryVector(n);
    int c = 0;
    for (auto &i : p.first) {
      z(c) = i;
      c++;
    }
    c = 0;
    for (auto &i : p.second) {
      x(c) = i;
      c++;
    }

    row << x, z;
    if (row != Eigen::VectorXi::Zero(2 * n)) {
      nRows++;
      tableaux.conservativeResize(nRows, 2 * n);
      tableaux.row(nRows - 1) = row;
    }
  }

  Eigen::MatrixXd tableauxD = tableaux.cast<double>();
  Eigen::MatrixXd b = tableauxD;
  std::vector<int> pivotCols;
  reduced_row_echelon_form(tableauxD, b, pivotCols);

  Eigen::MatrixXi linIndVecs(pivotCols.size(), 2 * n);
  linIndVecs.setZero();
  std::vector<std::vector<int>> linIndVecsVector;
  for (int i = 0; i < pivotCols.size(); i++) {
    linIndVecs.row(i) = b.cast<int>().row(i);
  }

  std::vector<int> zero(2 * n);
  std::set<std::vector<int>> generators, generated;
  generated.insert(zero);

  int counter = 0;
  std::vector<std::vector<int>> combinations;

  for (int nOnes = 0; nOnes <= n; nOnes++) {
    std::vector<int> test(n);

    for (int k = 0; k < nOnes; k++)
      test[n - k - 1] = 1;

    do {
      combinations.push_back(test);
      counter++;
    } while (std::next_permutation(test.begin(), test.end()));
  }

  auto bip = [](std::vector<int> &bv1, std::vector<int> &bv2) {
    std::vector<int> ax(bv1.begin(), bv1.begin() + bv1.size() / 2);
    std::vector<int> az(bv1.begin() + bv1.size() / 2, bv1.begin() + bv1.size());
    std::vector<int> bx(bv2.begin(), bv2.begin() + bv2.size() / 2);
    std::vector<int> bz(bv2.begin() + bv2.size() / 2, bv2.begin() + bv2.size());

    Eigen::VectorXi axv = Eigen::Map<Eigen::VectorXi>(ax.data(), ax.size());
    Eigen::VectorXi azv = Eigen::Map<Eigen::VectorXi>(az.data(), az.size());
    Eigen::VectorXi bxv = Eigen::Map<Eigen::VectorXi>(bx.data(), bx.size());
    Eigen::VectorXi bzv = Eigen::Map<Eigen::VectorXi>(bz.data(), bz.size());

    return (axv.dot(bzv) + azv.dot(bxv)) % 2;
  };

  counter = 0;
  int ker_dim = 2 * n - pivotCols.size();

  auto getUnion = [](const std::set<std::vector<int>> &a,
                     const std::set<std::vector<int>> &b) {
    std::set<std::vector<int>> result = a;
    result.insert(b.begin(), b.end());
    return result;
  };

  for (auto &gz : combinations) {
    std::vector<int> tau(n);
    tau.insert(tau.end(), gz.begin(), gz.end());

    if (generated.find(tau) != generated.end()) {
      continue;
    } else {

      int sum = 0;
      for (int i = 0; i < linIndVecs.rows(); i++) {
        Eigen::VectorXi row = linIndVecs.row(i);

        std::vector<int> h(2 * n);
        Eigen::VectorXi::Map(&h[0], 2 * n) = row;

        sum += bip(tau, h);
      }

      if (sum == 0) {
        generators.insert(tau);
        counter++;

        if (counter == ker_dim)
          break;

        for (auto &g : generators) {
          std::vector<int> s(2 * n);
          for (int i = 0; i < 2 * n; i++) {
            s[i] = (tau[i] + g[i]) % 2;
          }

          std::set<std::vector<int>> unionSet = getUnion(generators, generated);
          if (unionSet.find(s) == unionSet.end()) {
            generated.insert(s);
          }
        }
      }
    }
  }

  std::vector<PauliOperator> taus;
  for (auto &g : generators) {
    std::map<int, std::string> terms;
    for (int i = 0; i < 2 * n; i++) {
      if (g[i] == 1) {
        if (i < n) {
          terms.insert({i, "X"});
        } else {
          terms.insert({i - n, "Z"});
        }
      }
    }
    taus.emplace_back(terms);
  }

  auto had = [](const int i) {
    return (PauliOperator({{i, "X"}}) + PauliOperator({{i, "Z"}})) *
           (1.0 / std::sqrt(2));
  };

  auto hc = [](PauliOperator &op) {
    PauliOperator newOp = op;
    for (auto &kv : newOp) {
      kv.second.coeff() = std::conj(kv.second.coeff());
    }
    return newOp;
  };

  double invSqrt2 = 1.0 / std::sqrt(2.0);
  PauliOperator U(1.0);
  for (auto &t : taus) {
    int i = t.getTerms().begin()->second.ops().begin()->first;
    U *= (t + PauliOperator({{i, "X"}})) * invSqrt2;
    U *= had(i);
  }

  PauliOperator HPrime = hc(U) * H * U;

  nRows = 0;
  Eigen::MatrixXi tableaux2;
  for (auto &term : HPrime) {
    Eigen::VectorXi x(n), z(n), row(2 * n);
    auto p = term.second.toBinaryVector(n);
    int c = 0;
    for (auto &i : p.first) {
      z(c) = i;
      c++;
    }
    c = 0;
    for (auto &i : p.second) {
      x(c) = i;
      c++;
    }

    row << x, z;
    if (row != Eigen::VectorXi::Zero(2 * n)) {
      nRows++;
      tableaux2.conservativeResize(nRows, 2 * n);
      tableaux2.row(nRows - 1) = row;
    }
  }

  Eigen::MatrixXd tableauxD2 = tableaux2.cast<double>();
  Eigen::MatrixXd b2 = tableauxD2;
  std::vector<int> pivotCols2;
  reduced_row_echelon_form(tableauxD2, b2, pivotCols2);

  std::set<int> pivotSet(pivotCols2.begin(), pivotCols2.end()), doubleNSet,
      nSet;
  for (int i = 0; i < 2 * n; i++)
    doubleNSet.insert(i);
  for (int i = 0; i < n; i++)
    nSet.insert(i);

  std::set<int> phase_sites, keep_sites;
  std::set_difference(doubleNSet.begin(), doubleNSet.end(), pivotSet.begin(),
                      pivotSet.end(),
                      std::inserter(phase_sites, phase_sites.end()));

  std::set_difference(nSet.begin(), nSet.end(), phase_sites.begin(),
                      phase_sites.end(),
                      std::inserter(keep_sites, keep_sites.end()));

  counter = 0;
  std::vector<std::vector<int>> phase_configs;

  for (int nOnes = 0; nOnes <= phase_sites.size(); nOnes++) {
    std::vector<int> test(phase_sites.size());

    for (int k = 0; k < nOnes; k++)
      test[phase_sites.size() - k - 1] = 1;

    do {

      for (int i = 0; i < phase_sites.size(); i++) {
        if (test[i] == 0)
          test[i] = -1;
      }

      phase_configs.push_back(test);
      counter++;
    } while (std::next_permutation(test.begin(), test.end()));
  }

  PauliOperator actualReduced;
  double energy = 0.0;

  std::size_t dim = 1;
  std::size_t two = 2;
  for (int i = 0; i < n; i++)
    dim *= two;

  for (auto &phases : phase_configs) {
    counter = 0;
    std::map<int, int> Z2_dict;
    for (auto &z : phase_sites) {
      Z2_dict.insert({z, phases[counter]});
      counter++;
    }

    PauliOperator reduced;

    for (auto &kv : HPrime) {
      auto term = kv.second;
      if (term.ops().empty()) {
        reduced += PauliOperator(term.coeff());
      } else {
        std::map<int, std::string> newTerm;
        double phase = 1.0;

        for (auto &termKv : term.ops()) {
          if (keep_sites.find(termKv.first) != keep_sites.end()) {
            newTerm.insert({termKv.first, termKv.second});
          }

          if (Z2_dict.count(termKv.first)) {
            phase *= Z2_dict[termKv.first];
          }
        }
        reduced += PauliOperator(newTerm, phase * term.coeff());
      }
    }

    // if (reduced.getTerms().size() == 1 &&
    //     reduced.getTerms().begin()->second.ops().empty()) {
    //   continue;
    // }
    // get ground state energy of reduced

    auto getBitStrForIdx = [&](std::uint64_t i) {
      std::stringstream s;
      for (int k = n - 1; k >= 0; k--)
        s << ((i >> k) & 1);
      return s.str();
    };

    Eigen::MatrixXcd A(dim, dim);
    A.setZero();
    for (std::uint64_t myRow = 0; myRow < dim; myRow++) {
      auto rowBitStr = getBitStrForIdx(myRow);
      auto results = reduced.computeActionOnBra(rowBitStr);
      for (auto &result : results) {
        std::uint64_t k = std::stol(result.first, nullptr, 2);
        A(myRow, k) += result.second;
      }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
    auto reducedEnergy = es.eigenvalues()[0];
    // std::cout << "ENERGY: " << es.eigenvalues()[0] << "\n";

    if (reducedEnergy < energy) {
      energy = reducedEnergy;
      actualReduced = reduced;
    }
  }

  std::cout << "Actual Reduced: \n" << actualReduced.toString() << "\n";

  // Here I expect the ansatz to be functions[0]->getInstruction(0);
  auto ansatz = std::dynamic_pointer_cast<Function>(
      ir->getKernels()[0]->getInstruction(0));

  auto newIR = actualReduced.toXACCIR();

  if (ansatz) {
    for (auto &k : newIR->getKernels()) {
      k->insertInstruction(0, ansatz);
    }
  }

  return newIR;
}

} // namespace vqe
} // namespace xacc
