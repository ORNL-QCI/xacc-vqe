#include "RDMPurificationDecorator.hpp"
#include "FermionCompiler.hpp"
#include "IRProvider.hpp"
#include "InstructionIterator.hpp"
#include "PauliOperator.hpp"
#include "RDMGenerator.hpp"
#include "XACC.hpp"
#include <unsupported/Eigen/CXX11/TensorSymmetry>
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

  Eigen::array<int, 2> cc2({1, 3});
  Eigen::DynamicSGroup rho_pqrs_Sym;
  rho_pqrs_Sym.addAntiSymmetry(0, 1);
  rho_pqrs_Sym.addAntiSymmetry(2, 3);

  xacc::info("Filtering 2-RDM");
  Eigen::Tensor<std::complex<double>, 4> filtered_rhopqrs(nQubits, nQubits,
                                                          nQubits, nQubits);
  filtered_rhopqrs.setZero();
  for (int p = 0; p < nQubits; p++) {
    for (int q = 0; q < nQubits; q++) {
      for (int r = 0; r < nQubits; r++) {
        for (int s = 0; s < nQubits; s++) {
          if ((p <= q) && (r <= s)) {
            filtered_rhopqrs(p, q, r, s) = rho_pqrs(p, q, r, s);
          }
        }
      }
    }
  }

  //   PURIFY these...
  Eigen::Tensor<std::complex<double>, 4> rhopqrs_sq =
      filtered_rhopqrs * filtered_rhopqrs;
  Eigen::Tensor<std::complex<double>, 4> diff_pqrs =
      rhopqrs_sq - filtered_rhopqrs;
  Eigen::Tensor<std::complex<double>, 4> diffsq_pqrs = diff_pqrs * diff_pqrs;
  //
  // checking trace of 1-rdm
  Eigen::Tensor<std::complex<double>, 2> rhopq_tensor =
      filtered_rhopqrs.trace(cc2);
  Eigen::Tensor<std::complex<double>, 0> rhopq_trace_tensor =
      rhopq_tensor.trace();
  auto rhopq_trace = std::real(rhopq_trace_tensor(0));
  xacc::info("Trace of 1-rdm after filter: " + std::to_string(rhopq_trace));
  // end checking trace

  //   Eigen::Tensor<std::complex<double>, 0> begin_trace_tensor2 =
  //   diffsq_pqrs.trace(); auto tr_pqrs = std::real(begin_trace_tensor2(0));

  //   for (int p = 0; p < nQubits; p++){
  //       for (int q = 0; q < nQubits; q++){
  //           for (int r = 0; r < nQubits; r++){
  //               for (int s = 0; s < nQubits; s++){
  //                       xacc::info(std::to_string(p) + ", " +
  //                       std::to_string(q) + ", " + std::to_string(r) + ", " +
  //                       std::to_string(s) + ", " +
  //                       std::to_string(std::real(filtered_rhopqrs(p, q, r,
  //                       s))));
  //               }
  //           }
  //       }
  //   }

  auto tr_pqrs = 0.0;
  for (int p = 0; p < nQubits; p++) {
    for (int q = 0; q < nQubits; q++) {
      tr_pqrs += std::real(filtered_rhopqrs(p, q, p, q));
    }
  }
  xacc::info("value of sum(rho_pqpq) before purification: " +
             std::to_string(tr_pqrs));

  xacc::info("Purifying 2-rdm");

  // Purify rho_pqrs
  int count = 0;
  while (count < 10) {

    filtered_rhopqrs = (3 * rhopqrs_sq) - (2 * filtered_rhopqrs * rhopqrs_sq);
    //     filtered_rhopqrs = (3 * rhopqrs_sq) - (2 *
    //     filtered_rhopqrs.contract(rhopqrs_sq, rank4_dims));
    //     //rhopqrs_sq.contract(filtered_rhopqrs, rank4_dims));

    //     auto tmp_val2 = 0.0;
    //     for (int p = 0; p < nQubits; p++){
    //       for (int q = 0; q < nQubits; q++){
    //           tmp_val2 += std::real(filtered_rhopqrs(p, q, p, q));
    //         }
    //     }
    //     xacc::info("value of sum(rho_pqpq) during purification:" +
    //     std::to_string(tmp_val2));

    //     Eigen::Tensor<std::complex<double>, 0> tmp_trace_tensor =
    //     filtered_rhopqrs.trace(); auto tmp_trace =
    //     std::real(tmp_trace_tensor(0)); xacc::info("Tr(rhopqrs) during
    //     purification: " + std::to_string(tmp_trace));

    //     filtered_rhopqrs = (1.0 / tmp_val2) * filtered_rhopqrs;

    rhopqrs_sq = filtered_rhopqrs * filtered_rhopqrs;
    diff_pqrs = rhopqrs_sq - filtered_rhopqrs;
    diffsq_pqrs = diff_pqrs * diff_pqrs;

    Eigen::Tensor<std::complex<double>, 0> trace_tensor = diffsq_pqrs.trace();
    tr_pqrs = std::real(trace_tensor(0));

    xacc::info("Iter: " + std::to_string(count));
    xacc::info("diffsq_pqrs trace: " + std::to_string(tr_pqrs));
    count += 1;
  }

  xacc::info("Reconstructing rhopqrs");
  // reconstruct rhopqrs using symmetry rules
  for (int p = 0; p < nQubits; p++) {
    for (int q = 0; q < nQubits; q++) {
      for (int r = 0; r < nQubits; r++) {
        for (int s = 0; s < nQubits; s++) {
          rho_pqrs_Sym(filtered_rhopqrs, p, q, r, s) =
              filtered_rhopqrs(p, q, r, s);
          rho_pqrs_Sym(filtered_rhopqrs, r, s, p, q) =
              filtered_rhopqrs(p, q, r, s);
        }
      }
    }
  }

  auto filtered_trace = 0.0;
  for (int p = 0; p < nQubits; p++) {
    for (int q = 0; q < nQubits; q++) {
      filtered_trace += std::real(filtered_rhopqrs(p, q, p, q));
    }
  }

  xacc::info("value of sum(rhopqpq): " + std::to_string(filtered_trace));
  rhopq_tensor = filtered_rhopqrs.trace(cc2);
  rhopq_trace_tensor = rhopq_tensor.trace();
  rhopq_trace = std::real(rhopq_trace_tensor(0));
  xacc::info("Tr(rhopq): " + std::to_string(rhopq_trace));

  xacc::info("Computing purified energy");
  // Compute the energy
  xacc::info("Energy before calculation: " + std::to_string(energy));

  for (int p = 0; p < nQubits; p++) {
    for (int q = 0; q < nQubits; q++) {
      for (int r = 0; r < nQubits; r++) {
        for (int s = 0; s < nQubits; s++) {
          energy += 0.5 * std::real(hpqrs(p, q, r, s) *
                                    (filtered_rhopqrs(p, q, s, r) +
                                     filtered_rhopqrs(r, s, q, p)));
        }
      }
    }
  }
  xacc::info("Energy after 2-rdm: " + std::to_string(energy));

  for (int p = 0; p < nQubits; p++) {
    for (int q = 0; q < nQubits; q++) {
      energy += 0.5 * std::real(hpq(p, q) *
                                (rhopq_tensor(p, q) + rhopq_tensor(q, p)));
    }
  }
  xacc::info("Purified energy " + std::to_string(energy));

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
