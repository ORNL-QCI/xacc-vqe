#include "HWE.hpp"
#include "GateFunction.hpp"
#include <boost/tokenizer.hpp>

using namespace xacc::quantum;

namespace xacc {
namespace vqe {

std::shared_ptr<Function>
HWE::generate(std::map<std::string, InstructionParameter> parameters) {

  std::vector<InstructionParameter> params;
  if (!parameters.count("layers")) {
    params.push_back(InstructionParameter(1));
  }

  if (!parameters.count("n_qubits")) {
    xacc::error("HWE needs n_qubits parameter.");
  }

  params.push_back(parameters["n_qubits"]);

  if (!parameters.count("connectivity")) {
    xacc::error("HWE needs to know the connectivity of the qubits.");
  }

  params.push_back(parameters["connectivity"]);

  return generate(nullptr, params);
}

std::shared_ptr<Function>
HWE::generate(std::shared_ptr<AcceleratorBuffer> buffer,
              std::vector<InstructionParameter> parameters) {

  int nQubits = 0, layers = 1;
  std::vector<std::pair<int, int>> connectivity;
  try {
    nQubits = boost::get<int>(parameters[1]);
    layers = boost::get<int>(parameters[0]);
    auto cStr = boost::get<std::string>(parameters[2]);

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

  } catch (std::exception &e) {
    xacc::error("Could not cast HWE parameter to correct type: " +
                std::string(e.what()));
  }

  std::vector<InstructionParameter> fParams;
  for (int nP = 0; nP < 2 * nQubits * layers; nP++)
    fParams.push_back(InstructionParameter("t" + std::to_string(nP)));

  auto provider = xacc::getService<IRProvider>("gate");
  auto f = provider->createFunction("hwe", {}, fParams);
  int angleCounter = 0;
  for (int d = 0; d < layers; d++) {

    for (int q = 0; q < nQubits; q++) {
      auto rx = provider->createInstruction(
          "Rx", {q},
          {InstructionParameter("t" + std::to_string(angleCounter))});
      auto rz = provider->createInstruction(
          "Rz", {q},
          {InstructionParameter("t" + std::to_string(angleCounter + 1))});
      f->addInstruction(rx);
      f->addInstruction(rz);
      angleCounter += 2;
    }

    for (auto &p : connectivity) {
      auto cnot = provider->createInstruction("CNOT", {p.first, p.second});
      f->addInstruction(cnot);
    }
  }

  return f;
}

} // namespace vqe
} // namespace xacc
