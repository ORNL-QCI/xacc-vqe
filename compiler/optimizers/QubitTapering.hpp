#ifndef COMPILER_QUBIT_TAPERING_HPP_
#define COMPILER_QUBIT_TAPERING_HPP_

#include "GateFunction.hpp"
#include "IRTransformation.hpp"
#include "InstructionIterator.hpp"

namespace xacc {
namespace vqe {

class QubitTapering : public IRTransformation, public OptionsProvider {

public:
  QubitTapering() {}

  virtual std::shared_ptr<IR> transform(std::shared_ptr<IR> ir);

  virtual const std::string name() const { return "qubit-tapering"; }

  virtual const std::string description() const { return ""; }

  virtual std::shared_ptr<options_description> getOptions() {
    auto desc = std::make_shared<options_description>("QubitTapering Options");
    desc->add_options()("circuit-opt-n-tries", value<std::string>(),
                        "Provide the number of passes to use in optimizing "
                        "this circuit. Default = 2.")("circuit-opt-silent",
                                                      "Do not print any info");
    return desc;
  }

  virtual bool handleOptions(variables_map &map) { return false; }
};
} // namespace vqe
} // namespace xacc

#endif
