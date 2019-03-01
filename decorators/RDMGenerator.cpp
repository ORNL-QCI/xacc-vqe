#include "RDMGenerator.hpp"
#include <unsupported/Eigen/CXX11/TensorSymmetry>
#include "XACC.hpp"

namespace xacc {
namespace vqe {

void RDMGenerator::generate(std::shared_ptr<Function> ansatz) {
  // Reset
  rho_pq.setZero();
  rho_pqrs.setZero();

  // Get the Accelerator we're running on, the
  // number of orbitals/qubits in the problem,
  // the MPIProvider Communicator reference, and the
  // VQE state preparation ansatz.
  int nQubits = _nQubits, nExecs = 0;

  // Silence the FermionCompiler info messages.
  xacc::setOption("fermion-compiler-silent", "");

  // Get the FermionCompiler, we'll use it to map our
  // RDM second quantized expressions to Paulis that
  // can be measured on the QPU.
  auto fermionCompiler = xacc::getCompiler("fermion");

  Eigen::DynamicSGroup rho_pq_Sym, rho_pqrs_Sym;
  rho_pq_Sym.addHermiticity(0, 1);
  rho_pqrs_Sym.addAntiSymmetry(0, 1);
  rho_pqrs_Sym.addAntiSymmetry(2, 3);

  // Compute the 1 body reduced density matrix
  for (int m = 0; m < nQubits; m++) {
    for (int n = m; n < nQubits; n++) {
      xacc::info("Computing 1-rdm " + std::to_string(m) + ", " + std::to_string(n));
      // Generate the XACC kernel source code and compile it
      std::stringstream ss;
      ss << "__qpu__ k(){\n0.5 " << m << " 1 " << n << " 0\n"
         << "0.5 " << n << " 1 " << m << " 0\n}";
      double sum_pq = 0.0;
      auto hpq_ir = fermionCompiler->compile(ss.str(), qpu);

      // Loop over our compiled kernels and execute them
      // computing the weight*expVal sum
      for (auto &kernel : hpq_ir->getKernels()) {
        double localExpVal = 1.0;
        auto t = std::real(
            boost::get<std::complex<double>>(kernel->getParameter(0)));
        if (kernel->nInstructions() > 0) {

          kernel->insertInstruction(0, ansatz);

          auto buffer = qpu->createBuffer("q", nQubits);
          qpu->execute(buffer, kernel);
          localExpVal = buffer->getExpectationValueZ();
          nExecs++;
        }
        sum_pq += t * localExpVal;
      }

      if (std::fabs(sum_pq) > 1e-12)
        rho_pq_Sym(rho_pq, m, n) = sum_pq;
    }
  }

  double trace = 0.0;
  for (int i = 0; i < nQubits; i++)
    trace += std::real(rho_pq(i, i));

  // std::cout << "Trace for 1 RDM: " << trace << "\n";

  // Generate the 2-RDM
  for (int m = 0; m < nQubits; m++) {
    for (int n = m + 1; n < nQubits; n++) {
      for (int v = m; v < nQubits; v++) {
        for (int w = v + 1; w < nQubits; w++) {
      xacc::info("Computing 2-rdm " + std::to_string(m) + ", " + std::to_string(n) + ", " + std::to_string(v) +", " + std::to_string(w));

          // Create the source code XACC kernel and compile it
          std::stringstream xx;
          xx << "__qpu__ k(){\n0.5 " << m << " 1 " << n << " 1 " << w << " 0 "
             << v << " 0\n"
             << "0.5 " << w << " 1 " << v << " 1 " << m << " 0 " << n
             << " 0\n}";
          auto hpqrs_ir = fermionCompiler->compile(xx.str(), qpu);

          // Loop over the kernels, execute them and compute
          // the weight * expVal sum
          double sum_pqrs = 0.0;
          for (auto &kernel : hpqrs_ir->getKernels()) {
            double localExpVal = 1.0;
            auto t = std::real(
                boost::get<std::complex<double>>(kernel->getParameter(0)));
            if (kernel->nInstructions() > 0) {
              kernel->insertInstruction(0, ansatz);
              auto buffer = qpu->createBuffer("q", nQubits);
              qpu->execute(buffer, kernel);
              localExpVal = buffer->getExpectationValueZ();
              nExecs++;
            }

            sum_pqrs += t * localExpVal;
          }

          // Set the 2-RDM elements and update the energy.
          if (std::fabs(sum_pqrs) > 1e-12) {
            rho_pqrs_Sym(rho_pqrs, m, n, v, w) = sum_pqrs;
            rho_pqrs_Sym(rho_pqrs, v, w, m, n) = sum_pqrs;
          }
        }
      }
    }
  }

//   trace = 0;
//   for (int i = 0; i < nQubits; i++) {
//     for (int j = 0; j < nQubits; j++) {
//       trace += std::real(rho_pqrs(i, j, i, j));
//     }
//   }

//   // std::cout << "2 body trace = " << trace << "\n";
//   for (int m = 0; m < nQubits; m++) {
//     for (int n = 0; n < nQubits; n++) {
//       _energy += 0.5 * std::real(hpq(m, n) * (rho_pq(m, n) + rho_pq(n, m)));
//     }
//   }

//   for (int m = 0; m < nQubits; m++) {
//     for (int n = 0; n < nQubits; n++) {
//       for (int v = 0; v < nQubits; v++) {
//         for (int w = 0; w < nQubits; w++) {
//           _energy +=
//               0.5 * std::real(hpqrs(m, n, v, w) *
//                               (rho_pqrs(m, n, w, v) + rho_pqrs(v, w, n, m)));
//         }
//       }
//     }
//   }
  // std::cout << "RHOPQ:\n" << rho_pq << "\n";

  // for (int m = 0; m < nQubits; m++) {
  // 	for (int n = 0; n < nQubits; n++) {
  // 		for (int v = 0; v < nQubits; v++) {
  // 			for (int w = 0; w < nQubits; w++) {
  //                 if (std::fabs(rho_pqrs(m,n,v,w)) > 1e-12) std::cout <<
  //                 "2RDM: " << m << ", " << n << ", " << v << ", " << w << ":
  //                 " << rho_pqrs(m,n,v,w) << "\n";
  //             }
  //         }
  //     }
  // }

  // xacc::info("Number of QPU Executions = " + std::to_string(nExecs));
  // xacc::info("RDM generation complete.");

}

const double RDMGenerator::energy() { return _energy; }

} // namespace vqe
} // namespace xacc
