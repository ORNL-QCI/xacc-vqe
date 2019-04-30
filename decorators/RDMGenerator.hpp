#ifndef RDMGENERATOR_HPP_
#define RDMGENERATOR_HPP_

#include "Accelerator.hpp"
#include <unsupported/Eigen/CXX11/Tensor>

namespace xacc {
namespace vqe {

/**
 * The RDMGenerator provides a generate routine
 * that will create 1,2-RDMs for the fermionic Hamiltonian described
 * by the given VQEProgram.
 */
class RDMGenerator : public OptionsProvider {

protected:
  std::shared_ptr<Accelerator> qpu;

  /**
   * The energy computed by the RDM generation
   * routine.
   */
  double _energy = 0.0;

  int _nQubits;

public:
  /**
   * The 1-body RDM.
   */
  Eigen::Tensor<std::complex<double>, 2> rho_pq;

  /**
   * The 2-body RDM
   */
  Eigen::Tensor<std::complex<double>, 4> rho_pqrs;

  /**
   * The 1-body electron overlap integrals
   */
  Eigen::Tensor<std::complex<double>, 2> hpq;

  /**
   * The 2-body electron overlap integrals
   */
  Eigen::Tensor<std::complex<double>, 4> hpqrs;

  /**
   * The Constructor, takes a built VQEProgram
   * describing the high-level second quantized
   * fermionic hamiltonian.
   *
   * @param p
   */
  RDMGenerator(const int nQubits, std::shared_ptr<Accelerator> acc,
               Eigen::Tensor<std::complex<double>, 2> &hpq_,
               Eigen::Tensor<std::complex<double>, 4> &hpqrs_)
      : _nQubits(nQubits), rho_pq(nQubits, nQubits), qpu(acc),
        rho_pqrs(nQubits, nQubits, nQubits, nQubits), hpq(hpq_), hpqrs(hpqrs_) {
  }

  const int nQubits() { return _nQubits; }

  std::vector<std::shared_ptr<AcceleratorBuffer>> generate(std::shared_ptr<Function> ansatz, std::vector<int> qubitMap);

  /**
   * Return the computed hamiltonian energy.
   *
   * @return
   */
  const double energy();
    
  OptionPairs getOptions() override {
    OptionPairs desc {{"rdm-spins", ""}};
    return desc;
  }

};
} // namespace vqe
} // namespace xacc
#endif
