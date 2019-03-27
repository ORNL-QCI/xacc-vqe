#ifndef COMPILER_QUBIT_TAPERING_HPP_
#define COMPILER_QUBIT_TAPERING_HPP_

#include "IRTransformation.hpp"
#include "PauliOperator.hpp"
#include "OptionsProvider.hpp"
#include <Eigen/Dense>

using namespace xacc::quantum;

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

  virtual OptionPairs getOptions() {
    OptionPairs desc {{"phase-sector",
                        "Provide the +-1 vector."},{
        "qubit-tapering-show",
        "Create and display reduced hamiltonian, but then exit."}};
    return desc;
  }

private:
  Eigen::MatrixXi computeTableaux(PauliOperator &H, const int nQubits);
  std::vector<std::vector<int>> generateCombinations(
      const int nQubits, std::function<void(std::vector<int> &)> &&f =
                             [](std::vector<int> &tmp) { return; });
  int binaryVectorInnerProduct(std::vector<int> &bv1, std::vector<int> &bv2);
  const double computeGroundStateEnergy(PauliOperator &op, const int n);

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
};
} // namespace vqe
} // namespace xacc

#endif
