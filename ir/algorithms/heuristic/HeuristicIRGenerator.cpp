#include "HeuristicIRGenerator.hpp"
#include "GateFunction.hpp"
#include <boost/tokenizer.hpp>

using namespace xacc::quantum;

namespace xacc {
namespace vqe {

std::shared_ptr<Function>
HeuristicIRGenerator::generate(std::map<std::string, InstructionParameter> parameters) {

  std::vector<InstructionParameter> params;

  // Connectivity of the Hamiltonian
  if (!parameters.count("connectivity")) {
    xacc::error("HeuristicIRGenerator needs to know the connectivity of the Hamiltonian.");
  }

  params.push_back(parameters["connectivity"]);

  if (!parameters.count("interactions")) {
      xacc::error("HeuristicIRGenerator needs to know the interactions for each connection.");
  }
  
  params.push_back(parameters["interactions"]);
  
  return generate(nullptr, params);
}

std::shared_ptr<Function>
HeuristicIRGenerator::generate(std::shared_ptr<AcceleratorBuffer> buffer,
              std::vector<InstructionParameter> parameters) {

  std::vector<std::pair<int, int>> connectivity;
  std::vector<std::string> interactions;

  try {
    auto cStr = boost::get<std::string>(parameters[0]);
    auto iStr = boost::get<std::string>(parameters[1]);

    boost::replace_all(cStr, "'", "");
    boost::replace_all(iStr, "'", "");

    std::vector<int> ints;
    boost::char_separator<char> sep(",");
    boost::tokenizer<char_separator<char>> tokens(cStr, sep);
    for (auto &t : tokens) {
      std::stringstream s;
      s << t;
      auto tmp = s.str();
      boost::replace_all(tmp, "[", "");
      boost::replace_all(tmp, "]", "");
      boost::trim(tmp);
      int value;
      if (std::stringstream(tmp) >> value)
        ints.push_back(value);
    }

    if (ints.size() % 2 != 0)
      xacc::error(std::to_string(ints.size()) +
                  ", You must provide the connectivity as a list of pairs.");

    for (int i = 0; i < ints.size(); i += 2) {
      connectivity.push_back({ints[i], ints[i + 1]});
    }

    boost::tokenizer<char_separator<char>> interactionTokens(iStr, sep);
    for (auto &t : tokens) {
      std::stringstream s;
      s << t;
      auto tmp = s.str();
      boost::trim(tmp);
      interactions.push_back(tmp);
    }

  } catch (std::exception &e) {
    xacc::error("Could not cast HWE parameter to correct type: " +
                std::string(e.what()));
  }

  // You now have a list of interaction strings and 
  // a list of graph connections (connectivity and interactions vectors)
  auto provider = xacc::getService<IRProvider>("gate");

  // Create the Function to return. This needs to know the name 
  // of hte function, and a list of InstructionParameters that 
  // are the variable names for all the variational angles
  // Here's how its done in HWE
  // std::vector<InstructionParameter> fParams;
  // for (int nP = 0; nP < (nQubits + 3 * nQubits * layers); nP++)
        // fParams.push_back(InstructionParameter("t" + std::to_string(nP)));
  // auto f = provider->createFunction("hwe", {}, fParams);

 
  // Now just populate the Function with the addInstruction() call
  // Instructions are created with the provider->createInstruction() method
  // CREATE Rx example: auto rx = provider->createInstruction(
  //                         "Rx", {q}, {InstructionParameter("t" + std::to_string(angleCounter))});

  int angleCounter = 0;

  
  return nullptr;
}

} // namespace vqe
} // namespace xacc
