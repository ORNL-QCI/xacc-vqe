#include "VQERestartDecorator.hpp"
#include "IRProvider.hpp"
#include "PauliOperator.hpp"
#include "XACC.hpp"
#include <fstream>

namespace xacc {
namespace vqe {
void VQERestartDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<Function> function) {
  xacc::error("SymVerificationDecorator not implemented for single Function "
              "execution.");
  return;
}

void VQERestartDecorator::initialize() {

  if (!xacc::optionExists("vqe-restart-file")) {
    xacc::error("Cannot use VQERestartDecorator without vqe-restart-file option.");
  }
  auto fileStr = xacc::getOption("vqe-restart-file");
  std::ifstream t(fileStr);
  std::string restartABFile((std::istreambuf_iterator<char>(t)),
                   std::istreambuf_iterator<char>());

  std::istringstream is(restartABFile);
  loadedBuffer = std::make_shared<AcceleratorBuffer>();
  loadedBuffer->load(is);

  // So now we have a buffer that has N parameters,
  // and each set of parameters has M kernel children.
  // We should add each list of M into a queue, each time execute is
  // called we pop the list of children off the queue and return them
  // once the queue is empty, we just delegate execution

  // These should be given in the order they were produced
  auto params = loadedBuffer->getAllUnique("parameters");
  for (auto & p : params) {
      std::vector<std::shared_ptr<AcceleratorBuffer>> toAdd;
      auto children = loadedBuffer->getChildren("parameters",p);
      for (auto& c : children) {
          if ("I" != mpark::get<std::string>(c->getInformation("kernel"))) {
              toAdd.push_back(c);
          }
      }
      queuedChildren.push(toAdd);
  }
}

std::vector<std::shared_ptr<AcceleratorBuffer>>
VQERestartDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<Function>> functions) {

  std::vector<std::shared_ptr<AcceleratorBuffer>> buffers;

  // Must initialize here since we overrode (is that a word?)
  // initialize() and did not have a set decorated accelerator
//   if (decoratedAccelerator && !initialized) {
//       decoratedAccelerator->initialize();
//       initialized = true;
//   }

  if (!queuedChildren.empty()) {
      buffers = queuedChildren.front();
      queuedChildren.pop();
      xacc::info(std::to_string(queuedChildren.size() )+ ", [VQERestart] Returning Queued Child List of size " + std::to_string(buffers.size())+".");
      return buffers;
  } else {
      if (!initialized) decoratedAccelerator->initialize();

      return decoratedAccelerator->execute(buffer, functions);
  }
}

} // namespace vqe
} // namespace xacc