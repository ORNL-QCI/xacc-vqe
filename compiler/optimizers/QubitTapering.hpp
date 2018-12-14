#ifndef COMPILER_QUBIT_TAPERING_HPP_
#define COMPILER_QUBIT_TAPERING_HPP_

#include "IRTransformation.hpp"
#include "PauliOperator.hpp"

namespace xacc {
namespace vqe {

/**
 * QubitTapering is an IRTransformation that implements
 * the Hamiltonian reduction scheme for discrete
 * Z2 symmetries as described in https://arxiv.org/pdf/1701.08213.pdf.
 */
class QubitTapering : public IRTransformation, public OptionsProvider {

public:
  QubitTapering() {}

  virtual std::shared_ptr<IR> transform(std::shared_ptr<IR> ir);

  virtual const std::string name() const { return "qubit-tapering"; }
  virtual const std::string description() const {
    return "Reduce number of qubits required by exploiting Z2 symmetries.";
  }

  virtual bool handleOptions(variables_map &map) { return false; }
  virtual std::shared_ptr<options_description> getOptions() {
    auto desc = std::make_shared<options_description>("QubitTapering Options");
    desc->add_options()("phase-sector", value<std::string>(),
                        "Provide the +-1 vector.")(
        "qubit-tapering-show",
        "Create and display reduced hamiltonian, but then exit.");
    return desc;
  }

private:
  Eigen::MatrixXi computeTableaux(PauliOperator &H, const int nQubits);
  std::vector<std::vector<int>> generateCombinations(
      const int nQubits, std::function<void(std::vector<int> &)> &&f =
                             [](std::vector<int> &tmp) { return; });
  int binaryVectorInnerProduct(std::vector<int> &bv1, std::vector<int> &bv2);
  const double computeGroundStateEnergy(PauliOperator &op, const int n);

  unsigned pivot(Eigen::MatrixXd &B, unsigned r, unsigned c) {
    for (unsigned k = r; k < B.rows(); k++) {
      if (std::fabs(B(k, c)) > 1e-12) {
        return k;
      }
    }
    return B.rows();
  }
  
  Eigen::MatrixXi gauss(Eigen::MatrixXi& A, std::vector<int> &pivotCols) {
     int n = A.rows();
     int m = A.cols();

     int sc = 0;

     for (int i = 0; i < m; i++) {
         bool found_row = false;
         int ip = i-sc;
         for (int k = ip; k < n; k++) {
             if (A(k,i) == 1) {
                 found_row = true;
                 
                 if (k > ip) {
                     for (int j = i; j < m; j++) {
                         auto tmp = A(k,j);
                         A(k,j) = A(ip,j);
                         A(ip,j) = tmp;
                     }
                 }

                 for (int l = ip+1; l < n; l++) {
                     if (A(l,i) != 0) {
                         for (int j = ip; j < m; j++) {
                             A(l,j) = (A(l,j)+A(ip,j)) %2;
                         }
                     }
                 }
                 break;
             }
         }
         if (!found_row) sc += 1;
     }

    unsigned row2 = 0;
    for (unsigned col = 0; col < A.cols() && row2 < A.rows(); col++) {
      if (std::fabs(A(row2, col)) < 1e-12)
        continue;

      pivotCols.push_back(col);
      row2++;
    }
    
     return A;

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
    //   for (int m = 0; m < B.cols(); m++) B.row(index)(m) = (int) B.row(index)(m) % 2;

      for (j = 0; j < row; j++) {
        if (j == index)
          continue;

        B.row(j) += (-1.0 * B(j, i)) * B.row(index);
        // for (int m = 0; m < B.cols(); m++) B.row(j)(m) = (int) B.row(j)(m) % 2;
      }

      Eigen::MatrixXi tmp = B.cast<int>();
      B = B.cast<int>().unaryExpr([](const int x) { return x%2; }).cast<double>();
      
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
};
} // namespace vqe
} // namespace xacc

#endif
