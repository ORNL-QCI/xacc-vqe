#include "SymVerificationDecorator.hpp"
#include "IRProvider.hpp"
#include "PauliOperator.hpp"
#include "XACC.hpp"
#include "xacc_service.hpp"

using namespace xacc::quantum;

namespace xacc {
namespace vqe {
void SymVerificationDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<Function> function) {
  xacc::error("SymVerificationDecorator not implemented for single Function "
              "execution.");
  return;
}

std::vector<std::shared_ptr<AcceleratorBuffer>>
SymVerificationDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<Function>> functions) {

  std::vector<std::shared_ptr<AcceleratorBuffer>> buffers;
  std::vector<std::shared_ptr<Function>> notConstFunctions;

  for (auto &f : functions)
    notConstFunctions.push_back(f);

  if (!xacc::optionExists("sym-op")) {
    xacc::error("Cannot use SymVerificationDecorator without sym-op option.");
  }

  int s = -1;
  if (xacc::optionExists("sym-s")) {
    s = std::stoi(xacc::getOption("sym-s"));
  }

  PauliOperator S;
  S.fromString(xacc::getOption("sym-op"));
  auto SName = S.getTerms().begin()->first;

  // Extract common instructions from the vector of Functions,
  // this makes up our ansatz
  auto ansatz =
      std::dynamic_pointer_cast<Function>(functions[0]->getInstruction(0));

  PauliOperator H;
  auto ir = xacc::getService<xacc::IRProvider>("gate")->createIR();
  for (auto &f : functions) {
    ir->addKernel(f);
  }
  H.fromXACCIR(ir);

  // Ensure that S commutes with H
  if (!H.commutes(S)) {
    xacc::error("The provided S is not a symmetry of H (S = " + S.toString() +
                ").");
  }

  // See if we need to add S to the functions list
  if (!H.contains(S)) {
    xacc::info("Adding S to Measure Kernels.");
    auto SFunction = S.toXACCIR()->getKernels()[0];
    SFunction->insertInstruction(0, ansatz);
    notConstFunctions.push_back(SFunction);
  }

  // See if we need to add PS to the functions list
  // for each P in H
  std::map<std::string, std::string> pToPS;
  std::map<std::string, double> coeffMap;
  for (auto &term : H) {
    PauliOperator P(term.second.ops());

    auto PS = P * S;
    auto tmpName = PS.getTerms().begin()->first;
    pToPS.insert({term.first, tmpName});
    coeffMap.insert(
        {tmpName, std::real(PS.getTerms().begin()->second.coeff())});

    if (!H.contains(PS) && PS.nTerms() > 0) {
      auto PSFunction = PS.toXACCIR()->getKernels()[0];
      if (PSFunction->nInstructions() > 0) {
        xacc::info("P * S\n" + P.toString() + " * " + S.toString() + " = " +
                   PS.toString());
        PSFunction->insertInstruction(0, ansatz);
        notConstFunctions.push_back(PSFunction);
      }
    }
  }

  // Execute with notConstFunctions
  xacc::info(notConstFunctions[0]->toString("q"));
  auto tmpBuffers = decoratedAccelerator->execute(buffer, notConstFunctions);

  bool useROEMExps = xacc::optionExists("sym-use-ro-error");

  // Create buffer name to expectation value map
  std::map<std::string, double> bufferMap;
  for (auto &b : tmpBuffers) {
    auto expval = (useROEMExps ? mpark::get<double>(b->getInformation("fixed-ro-error-exp-z")) : b->getExpectationValueZ());
    bufferMap.insert({b->name(), expval});
    xacc::info("BufMap: " + b->name() + " = " + std::to_string(expval));
  }

  double expS = bufferMap[SName];

  for (int i = 0; i < functions.size(); i++) {
    auto b = tmpBuffers[i];
    auto PName = b->name();
    auto PSName = pToPS[PName];

    double expPS = 1.0;
    if ("I" != PSName) {
      expPS = bufferMap[PSName];
    }

    auto expP = bufferMap[PName];

    auto fixed =
        (expP + s * expPS) / (1.0 + s * expS);

    xacc::info(PName + ", " + std::to_string(expP) + ", " +
               std::to_string(fixed) + ", " + std::to_string(expPS) + ", " +
               std::to_string(expS));

    b->addExtraInfo("sym-verification-fixed-exp-z", ExtraInfo(fixed));

    buffers.push_back(b);
  }

  buffer->addExtraInfo("sym-op", xacc::getOption("sym-op"));

  return buffers;
}

} // namespace vqe
} // namespace xacc