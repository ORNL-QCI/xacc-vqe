#include "RDMGenerator.hpp"
#include "XACC.hpp"
#include <unsupported/Eigen/CXX11/TensorSymmetry>

namespace xacc {
namespace vqe {

std::vector<std::shared_ptr<AcceleratorBuffer>> RDMGenerator::generate(std::shared_ptr<Function> ansatz, std::vector<int> qubitMap) {
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

  // Map p,q,r,s indices to a coefficient for all
  // identity terms encountered
  std::map<std::vector<int>, double> rho_element_2_identity_coeff;

  // Map function names to the actual function and all
  // p,q,r,s and coefficients that contribute to that element
  std::map<std::string,
           std::pair<std::shared_ptr<Function>,
                     std::vector<std::pair<std::vector<int>, double>>>>
      functions;

  // Generate the 2-RDM circuits for executiong
  for (int m = 0; m < nQubits; m++) {
    for (int n = m + 1; n < nQubits; n++) {
      for (int v = m; v < nQubits; v++) {
        for (int w = v + 1; w < nQubits; w++) {
          // Create the source code XACC kernel and compile it
          std::stringstream xx;
          xx << "__qpu__ k(){\n0.5 " << qubitMap[m] << " 1 " << qubitMap[n] << " 1 " << qubitMap[w] << " 0 "
             << qubitMap[v] << " 0\n"
             << "0.5 " << qubitMap[w] << " 1 " << qubitMap[v] << " 1 " << qubitMap[m] << " 0 " << qubitMap[n]
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
            //   xacc::info("kernel " + kernel->name() + ":\n" + kernel->toString());
              kernel->insertInstruction(0, ansatz);
              auto name = kernel->name();
              if (functions.count(name)) {
                functions[name].second.push_back({{m, n, v, w}, t});
              } else {
                functions.insert({name, {kernel, {{{m, n, v, w}, t}}}});
              }

            } else {
              rho_element_2_identity_coeff.insert({{m, n, v, w}, t});
            }
          }
        }
      }
    }
  }

  std::vector<std::shared_ptr<Function>> fsToExecute;
  for (auto &kv : functions) {
    fsToExecute.push_back(kv.second.first);
  }

  int nPhysicalQubits = *std::max_element(qubitMap.begin(), qubitMap.end()) + 1;

  // Execute all nontrivial circuits
  xacc::info(std::to_string(nPhysicalQubits) +", Executing " + std::to_string(fsToExecute.size()) +
             " circuits to compute rho_pqrs.");
  auto buffer = qpu->createBuffer("q", nPhysicalQubits);
  auto buffers = qpu->execute(buffer, fsToExecute);

  // Create a mapping of rho_pqrs elements to summed expectation values for
  // each circuit contributing to it
  std::map<std::vector<int>, double> sumMap;
  for (int i = 0; i < buffers.size(); i++) {
    auto fName = fsToExecute[i]->name();
    auto p = functions[fName];
    std::vector<std::string> contributingIndices;
    std::vector<double> contributingCoeffs;
    for (auto &l : p.second) {
      auto elements = l.first;
      std::stringstream s;
      s << elements[0] << "," << elements[1] << "," << elements[2] << "," << elements[3];
      contributingIndices.push_back(s.str());
      auto value = l.second * buffers[i]->getExpectationValueZ();
      contributingCoeffs.push_back(value);
      if (!sumMap.count(elements)) {
        sumMap.insert({elements, value});
      } else {
        sumMap[elements] += value;
      }
    }
    buffers[i]->addExtraInfo("kernel", ExtraInfo(fName));
    buffers[i]->addExtraInfo("contributing_rho_pqrs", ExtraInfo(contributingIndices));
    buffers[i]->addExtraInfo("contributing_coeffs", ExtraInfo(contributingCoeffs));
  }

  // Add all identity terms, we don't execute them
  // but we still have to add their contribution
  for (auto &kv : rho_element_2_identity_coeff) {
    sumMap[kv.first] += kv.second;
  }

  // Set rho_pqrs. This is all we need
  // to get rho_pq as well
  for (auto &kv : sumMap) {
    auto elements = kv.first;
    rho_pqrs_Sym(rho_pqrs, elements[0], elements[1], elements[2], elements[3]) =
        kv.second;
    rho_pqrs_Sym(rho_pqrs, elements[2], elements[3], elements[0], elements[1]) =
        kv.second;
  }

  return buffers;
}

const double RDMGenerator::energy() { return _energy; }

} // namespace vqe
} // namespace xacc
