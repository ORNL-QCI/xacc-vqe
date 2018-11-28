#include "PurificationDecorator.hpp"
#include "IRProvider.hpp"
#include "InstructionIterator.hpp"
#include "PauliOperator.hpp"
#include "XACC.hpp"

namespace xacc {
namespace quantum {
void PurificationDecorator::execute(std::shared_ptr<AcceleratorBuffer> buffer,
                                    const std::shared_ptr<Function> function) {

  if (!decoratedAccelerator) {
    xacc::error("Null Decorated Accelerator Error");
  }

  return;
}

std::vector<std::shared_ptr<AcceleratorBuffer>> PurificationDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<Function>> functions) {

  std::vector<std::shared_ptr<AcceleratorBuffer>> buffers;

  if (!decoratedAccelerator) {
    xacc::error("PurificationDecorator - Null Decorated Accelerator Error");
  }

  // Here I expect the ansatz to be functions[0]->getInstruction(0);
  auto ansatz =
      std::dynamic_pointer_cast<Function>(functions[0]->getInstruction(0));
      
  if (!ansatz)
    xacc::error("ANSATZ IS NULL");

  // Generate all nQubit Pauli Strings
  std::vector<std::string> XYZ;
  std::set<std::string> temp, PauliStrings;
  std::vector<xacc::vqe::PauliOperator> Paulis;

  for (int i = 0; i < buffer->size(); i++)
    XYZ.push_back("X" + std::to_string(i));
  for (int i = 0; i < buffer->size(); i++)
    XYZ.push_back("Y" + std::to_string(i));
  for (int i = 0; i < buffer->size(); i++)
    XYZ.push_back("Z" + std::to_string(i));

  std::function<void(std::vector<std::string> &, std::string, const int,
                     const int)>
      y;
  y = [&](std::vector<std::string> &set, std::string prefix, const int n,
          const int k) {
    if (k == 0) {
      xacc::vqe::PauliOperator op;
      op.fromString(prefix);
      PauliStrings.insert(op.toString());
      return;
    }

    for (int i = 0; i < n; i++) {
      auto newPrefix = prefix + " " + set[i];
      y(set, newPrefix, n, k - 1);
    }
  };

  auto x = [&](std::vector<std::string> &set, const int k) {
    auto n = set.size();
    y(set, std::string(""), n, k);
  };

  // Get all Permutations of XYZ
  x(XYZ, buffer->size());

  xacc::vqe::PauliOperator all;
  std::map<std::string, xacc::vqe::PauliOperator> opMap;
  for (auto &s : PauliStrings) {
    xacc::vqe::PauliOperator op;
    op.fromString(s);
    Paulis.push_back(op);
    all += op;
    opMap.insert({op.getTerms().begin()->second.id(), op});
  }

  for (auto &kv : opMap)
    std::cout << kv.first << ", " << kv.second.toString() << "\n";
  std::cout << all.toString() << "\n";
  std::cout << Paulis.size() << " operators to measure\n";

  auto kernels = all.toXACCIR()->getKernels();
  for (auto &k : kernels) {
    k->insertInstruction(0, ansatz);
  }

  std::cout << "EXECUTING\n";
  buffers = decoratedAccelerator->execute(buffer, kernels);

  std::size_t dim = 1;
  std::size_t two = 2;
  for (int i = 0; i < buffer->size(); i++)
    dim *= two;

  Eigen::MatrixXcd rho(dim, dim);
  rho.setZero();

  for (auto &b : buffers) {
    auto id = b->name();
    auto op = opMap[id];

    Eigen::MatrixXcd p(dim, dim);
    p.setZero();

    for (auto &t : op.getSparseMatrixElements()) {
      p(t.row(), t.col()) = t.coeff();
    }

    rho += b->getExpectationValueZ() * p;
  }

  // Add the Identity term
  rho += Eigen::MatrixXcd::Identity(dim, dim);

  rho *= 0.25;

  std::cout << "RHO:\n" << rho << "\n";
  std::cout << rho.trace() << "\n";

  Eigen::MatrixXcd rhosq = rho * rho;
  Eigen::MatrixXcd diff = rhosq - rho;
  Eigen::MatrixXcd diffsq = diff * diff;

  int counter = 0;
  auto tr = std::real(diffsq.trace());
  while (tr > 1e-2) {
    rho = 3. * rhosq - 2. * rhosq * rho;
    rho = rho / rho.trace();

    rhosq = rho * rho;
    diff = rhosq - rho;
    diffsq = diff * diff;
    tr = std::real(diffsq.trace());
    std::cout << counter << ", TRACE: " << tr << "\n";
    counter++;
    if (counter > 100)
      break;
  }

  std::cout << rho.trace() << "\n" << rho << "\n";

  // new E = Tr(H*rho)
  // so get H
  xacc::vqe::PauliOperator HOp;
  auto ir = xacc::getService<xacc::IRProvider>("gate")->createIR();
  for (auto &f : functions) {
    ir->addKernel(f);
  }
  HOp.fromXACCIR(ir);

  std::cout << "MADE IT HERE\n";
  std::cout << HOp.toString() << "\n";

  auto getBitStrForIdx = [&](std::uint64_t i) {
    std::stringstream s;
    for (int k = buffer->size() - 1; k >= 0; k--)
      s << ((i >> k) & 1);
    return s.str();
  };
  
  Eigen::MatrixXcd H(dim, dim);
  H.setZero();
  for (std::uint64_t myRow = 0; myRow < dim; myRow++) {
    auto rowBitStr = getBitStrForIdx(myRow);
    auto results = HOp.computeActionOnBra(rowBitStr);
    for (auto &result : results) {
      std::uint64_t k = std::stol(result.first, nullptr, 2);
      H(myRow, k) += result.second;
    }
  }

  H(0,0) = std::complex<double>(0.0,0.0);

  std::cout << "CONJ:\n" << (rho - rho.conjugate()) << "\n";
  std::cout << "IDEMP:\n" << (rho*rho - rho) << "\n";
  std::cout << "H:\n" << H << "\n";
  std::cout << "HI\n";
  std::cout << "Energy: " << (H * rho).trace() << "\n";
  return buffers;
}

} // namespace quantum
} // namespace xacc

// auto getBitStrForIdx = [&](std::uint64_t i) {
//     std::stringstream s;
//     for (int k = buffer->size() - 1; k >= 0; k--)
//       s << ((i >> k) & 1);
//     return s.str();
//   };
// for (std::uint64_t myRow = 0; myRow < dim; myRow++) {
//   auto rowBitStr = getBitStrForIdx(myRow);
//   auto results = op.computeActionOnBra(rowBitStr);
//   for (auto &result : results) {
//     std::uint64_t k = std::stol(result.first, nullptr, 2);
//     p(myRow, k) = result.second;
//   }
// }