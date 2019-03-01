#include "RDMPurificationDecorator.hpp"
#include "FermionCompiler.hpp"
#include "IRProvider.hpp"
#include "InstructionIterator.hpp"
#include "PauliOperator.hpp"
#include "RDMGenerator.hpp"
#include "XACC.hpp"
namespace xacc {
namespace vqe {
void RDMPurificationDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<Function> function) {

  if (!decoratedAccelerator) {
    xacc::error("Null Decorated Accelerator Error");
  }

  return;
}

std::vector<std::shared_ptr<AcceleratorBuffer>>
RDMPurificationDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<Function>> functions) {

  std::vector<std::shared_ptr<AcceleratorBuffer>> buffers;

  if (!decoratedAccelerator) {
    xacc::error("RDMPurificationDecorator - Null Decorated Accelerator Error");
  }

  // Here I expect the ansatz to be functions[0]->getInstruction(0);
  auto ansatz =
      std::dynamic_pointer_cast<Function>(functions[0]->getInstruction(0));

  if (!ansatz)
    xacc::error("ANSATZ IS NULL");

  auto src = xacc::getOption("rdm-source");

  auto nQubits = buffer->size();
  // Get hpq, hpqrs
  xacc::setOption("no-fermion-transformation", "");
  auto c = xacc::getCompiler("fermion");
  auto ir = c->compile(src, decoratedAccelerator);
  auto fermionKernel =
      std::dynamic_pointer_cast<FermionKernel>(ir->getKernels()[0]);
  xacc::unsetOption("no-fermion-transformation");

  auto energy = fermionKernel->E_nuc();
  auto hpq = fermionKernel->hpq(nQubits);
  auto hpqrs = fermionKernel->hpqrs(nQubits);
  RDMGenerator generator(buffer->size(), decoratedAccelerator, hpq, hpqrs);

  generator.generate(ansatz);

  auto rho_pq = generator.rho_pq;
  auto rho_pqrs = generator.rho_pqrs;

  // PURIFY these...
  Eigen::Tensor<std::complex<double>, 2> rhosq = rho_pq * rho_pq;
  Eigen::Tensor<std::complex<double>, 4> rhopqrs_sq = rho_pqrs * rho_pqrs;
  Eigen::Tensor<std::complex<double>, 2> diff = rhosq - rho_pq;
  Eigen::Tensor<std::complex<double>, 4> diff_pqrs = rhopqrs_sq - rho_pqrs;
  Eigen::Tensor<std::complex<double>, 2> diffsq = diff * diff;
  Eigen::Tensor<std::complex<double>, 4> diffsq_pqrs = diff_pqrs * diff_pqrs;

  double tr_pq = 0.0;
  for (int i = 0; i < nQubits; i++) {
    tr_pq += std::real(diffsq(i, i));
  }

  double tr_pqrs = 0;
  for (int i = 0; i < nQubits; i++) {
    for (int j = 0; j < nQubits; j++) {
      tr_pqrs += std::real(diffsq_pqrs(i, j, i, j));
    }
  }

  xacc::info("Purifying 1-rdm");
  std::cout << "HI: " << rho_pq << "\n";
  // Purify rho_pq
  while (tr_pq > 1e-4) {
    rho_pq = 3. * rhosq - 2. * rhosq * rho_pq;

    double tmp_trace = 0.;
    for (int i = 0; i < nQubits; i++) {
      tmp_trace += std::real(rho_pq(i, i));
    }
    rho_pq = (1.0 / tmp_trace) * rho_pq;

    rhosq = rho_pq * rho_pq;
    diff = rhosq - rho_pq;
    diffsq = diff * diff;
    tr_pq = 0.0;
    for (int i = 0; i < nQubits; i++) {
      tr_pq += std::real(diffsq(i, i));
    }
  }

    xacc::info("Purifying 2-rdm");


  // Purify rho_pqrs
  while (tr_pqrs > 1e-4) {
    rho_pqrs = 3. * rhopqrs_sq - 2. * rhopqrs_sq * rho_pqrs;

    double tmp_trace = 0.;
    for (int i = 0; i < nQubits; i++) {
      for (int j = 0; j < nQubits; j++) {
        tmp_trace += std::real(rho_pqrs(i, j, i, j));
      }
    }
    rho_pqrs = (1.0 / tmp_trace) * rho_pqrs;

    rhopqrs_sq = rho_pqrs * rho_pqrs;
    diff_pqrs = rhopqrs_sq - rho_pqrs;
    diffsq_pqrs = diff_pqrs * diff_pqrs;

    tr_pqrs = 0;
    for (int i = 0; i < nQubits; i++) {
      for (int j = 0; j < nQubits; j++) {
        tr_pqrs += std::real(diffsq_pqrs(i, j, i, j));
      }
    }
  }

  xacc::info("Computing purified energy");
  // Compute the energy
  for (int m = 0; m < nQubits; m++) {
    for (int n = 0; n < nQubits; n++) {
      energy += 0.5 * std::real(hpq(m, n) * (rho_pq(m, n) + rho_pq(n, m)));
    }
  }

  for (int m = 0; m < nQubits; m++) {
    for (int n = 0; n < nQubits; n++) {
      for (int v = 0; v < nQubits; v++) {
        for (int w = 0; w < nQubits; w++) {
          energy += 0.5 * std::real(hpqrs(m, n, v, w) * (rho_pqrs(m, n, w, v) +
                                                         rho_pqrs(v, w, n, m)));
        }
      }
    }
  }

  xacc::info("purified energy " + std::to_string(energy));

  std::vector<std::shared_ptr<AcceleratorBuffer>> retBuffers;
  for (auto &f : functions) {
    auto b = decoratedAccelerator->createBuffer(f->name(), buffer->size());
    b->addExtraInfo("purified-energy", ExtraInfo(energy));
    retBuffers.push_back(b);
  }

  return retBuffers;
}

} // namespace vqe
} // namespace xacc
