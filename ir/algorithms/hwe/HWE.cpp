#include "HWE.hpp"
#include "GateFunction.hpp"
// #include <boost/tokenizer.hpp>
#include "xacc_service.hpp"
#include <regex>

using namespace xacc::quantum;

namespace xacc {
namespace vqe {

std::shared_ptr<Function>
HWE::generate(std::map<std::string, InstructionParameter> parameters) {

  std::vector<InstructionParameter> params;
  if (!parameters.count("layers")) {
    params.push_back(InstructionParameter(1));
  } else {
    params.push_back(parameters["layers"]);
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
    nQubits = parameters[1].as<int>();
    layers = parameters[0].as<int>();
    auto cStr = parameters[2].as<std::string>();
    cStr = std::regex_replace(cStr, std::regex("'"), "");

    // boost::replace_all(cStr, "'", "");

    // xacc::info("CONNECTIVITY STRING: " + cStr);
    std::vector<int> ints;
    // boost::char_separator<char> sep(",");
    // boost::tokenizer<char_separator<char>> tokens(cStr, sep);
    auto tokens = xacc::split(cStr, ',');
    for (auto &t : tokens) {
      std::stringstream s;
      s << t;
      auto tmp = s.str();
      tmp = std::regex_replace(tmp, std::regex("\\["), "");
      tmp = std::regex_replace(tmp, std::regex("\\]"), "");

      xacc::trim(tmp);
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
  for (int nP = 0; nP < (2*nQubits + 3 * nQubits * layers); nP++)
    fParams.push_back(InstructionParameter("t" + std::to_string(nP)));

  auto provider = xacc::getService<IRProvider>("gate");
  auto f = provider->createFunction("hwe", {}, fParams);
  int angleCounter = 0;

  // Zeroth layer, start with X and Z rotations
  for (int q = 0; q < nQubits; q++) {
    auto rx = provider->createInstruction(
        "Rx", {q}, {InstructionParameter("t" + std::to_string(angleCounter))});
    auto rz = provider->createInstruction(
        "Rz", {q}, {InstructionParameter("t" + std::to_string(angleCounter+1))});
    f->addInstruction(rx);
    f->addInstruction(rz);
    angleCounter+=2;
  }

  for (int d = 0; d < layers; d++) {
    for (auto &p : connectivity) {
      auto cnot = provider->createInstruction("CNOT", {p.first, p.second});
      f->addInstruction(cnot);
    }
    for (int q = 0; q < nQubits; q++) {
      auto rz1 = provider->createInstruction(
          "Rz", {q},
          {InstructionParameter("t" + std::to_string(angleCounter))});
      f->addInstruction(rz1);

      auto rx = provider->createInstruction(
          "Rx", {q},
          {InstructionParameter("t" + std::to_string(angleCounter+1))});
      f->addInstruction(rx);

      auto rz2 = provider->createInstruction(
          "Rz", {q},
          {InstructionParameter("t" + std::to_string(angleCounter + 2))});
      f->addInstruction(rz2);

      angleCounter += 3;
    }
  }

  //   std::cout << "THE FUNCTION:\n" << f->toString("q") << "\n";
  return f;
}

} // namespace vqe
} // namespace xacc
