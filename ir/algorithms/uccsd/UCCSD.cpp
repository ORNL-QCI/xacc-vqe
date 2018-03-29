#include "UCCSD.hpp"
#include "GateFunction.hpp"
#include "FermionToSpinTransformation.hpp"
#include "ServiceRegistry.hpp"
#include "CommutingSetGenerator.hpp"
#include <boost/math/constants/constants.hpp>

using namespace xacc::quantum;

namespace xacc {
namespace vqe {


std::shared_ptr<Function> UCCSD::generate(
		std::shared_ptr<AcceleratorBuffer> buffer,
		std::vector<InstructionParameter> parameters) {

	auto runtimeOptions = RuntimeOptions::instance();

	if (!runtimeOptions->exists("n-electrons")) {
		xacc::error("To use this UCCSD State Prep IRGenerator, you "
				"must specify the number of electrons.");
	}

	if (!runtimeOptions->exists("n-qubits")) {
		xacc::error("To use this UCCSD State Prep IRGenerator, you "
				"must specify the number of qubits.");
	}

	auto nQubits = std::stoi((*runtimeOptions)["n-qubits"]);
	auto nElectrons = std::stoi((*runtimeOptions)["n-electrons"]);

	// Compute the number of parameters
	auto _nOccupied = (int) std::ceil(nElectrons / 2.0);
	auto _nVirtual = nQubits / 2 - _nOccupied;
	auto nSingle = _nOccupied * _nVirtual;
	auto nDouble = std::pow(nSingle, 2);
	auto _nParameters = nSingle + nDouble;

	auto singletIndex = [=](int i, int j) -> int {
		return i * _nOccupied + j;
	};

	auto doubletIndex = [=](int i, int j, int k, int l) -> int {
		return
		(i * _nOccupied * _nVirtual * _nOccupied +
				j * _nVirtual * _nOccupied +
				k * _nOccupied +
				l);
	};

	std::vector<xacc::InstructionParameter> variables;
	std::vector<std::string> params;
	for (int i = 0; i < _nParameters; i++) {
		params.push_back("theta" + std::to_string(i));
		variables.push_back(InstructionParameter("theta" + std::to_string(i)));
	}

	auto kernel = std::make_shared<FermionKernel>("fermiUCCSD");
	xacc::info("Constructing UCCSD Fermion Operator.");
	for (int i = 0; i < _nVirtual; i++) {
		for (int j = 0; j < _nOccupied; j++) {
			for (int l = 0; l < 2; l++) {
				std::vector<std::pair<int, int>> operators { { 2
						* (i + _nOccupied) + l, 1 }, { 2 * j + l, 0 } };
				auto fermiInstruction1 = std::make_shared<FermionInstruction>(
						operators, params[singletIndex(i, j)]);
				kernel->addInstruction(fermiInstruction1);

				std::vector<std::pair<int, int>> operators2 { { 2 * j + l, 1 },
						{ 2 * (i + _nOccupied) + l, 0 } };
				auto fermiInstruction2 = std::make_shared<FermionInstruction>(
						operators2, params[singletIndex(i, j)]);


				auto nP = fermiInstruction2->nParameters();
				InstructionParameter p(-1.0*std::complex<double>(1,0));
				fermiInstruction2->setParameter(nP-2,p);

				kernel->addInstruction(fermiInstruction2);

			}
		}
	}

	for (int i = 0; i < _nVirtual; i++) {
		for (int j = 0; j < _nOccupied; j++) {
			for (int l = 0; l < 2; l++) {
				for (int i2 = 0; i2 < _nVirtual; i2++) {
					for (int j2 = 0; j2 < _nOccupied; j2++) {
						for (int l2 = 0; l2 < 2; l2++) {
							std::vector<std::pair<int, int>> operators1 { { 2
									* (i + _nOccupied) + l, 1 },
									{ 2 * j + l, 0 }, { 2 * (i2 + _nOccupied)
											+ l2, 1 }, { 2 * j2 + l2, 0 } };

							std::vector<std::pair<int, int>> operators2 { { 2
									* j2 + l2, 1 }, { 2 * (i2 + _nOccupied)
									+ l2, 0 }, { 2 * j + l, 1 }, { 2
									* (i + _nOccupied) + l, 0 } };

							auto doubletIdx1 = nSingle
									+ doubletIndex(i, j, i2, j2);
							auto doubletIdx2 = nSingle
									+ doubletIndex(i, j, i2, j2);

							auto fermiInstruction1 = std::make_shared<
									FermionInstruction>(operators1,
									params[doubletIdx1]);

							kernel->addInstruction(fermiInstruction1);

							auto fermiInstruction2 = std::make_shared<
									FermionInstruction>(operators2,
									params[doubletIdx2]);

							auto nP = fermiInstruction2->nParameters();
							InstructionParameter p(-1.0*std::complex<double>(1,0));
							fermiInstruction2->setParameter(nP-2,p);
							kernel->addInstruction(fermiInstruction2);

						}
					}
				}
			}
		}
	}

//	std::cout << "KERNEL: \n" << kernel->toString("") << "\n";
	// Create the FermionIR to pass to our transformation.
	auto fermionir = std::make_shared<FermionIR>();
	fermionir->addKernel(kernel);

	xacc::info("Done constructing UCCSD Fermion Operator.");
	xacc::info("Mapping UCCSD Fermion Operator to Spin. ");

	std::shared_ptr<FermionToSpinTransformation> transform;
	if (xacc::optionExists("fermion-transformation")) {
		auto transformStr = xacc::getOption("fermion-transformation");
		transform = ServiceRegistry::instance()->getService<FermionToSpinTransformation>(
				transformStr);
	} else {
		transform = ServiceRegistry::instance()->getService<FermionToSpinTransformation>(
				"jw");
	}

	auto compositeResult = transform->transform(*kernel.get());
//	auto resultsStr = compositeResult.toString();
//	boost::replace_all(resultsStr, "+", "+\n");

	// Create the Spin Hamiltonian
	auto transformedIR = compositeResult.toXACCIR();
	xacc::info("Done mapping UCCSD Fermion Operator to Spin.");

	std::unordered_map<std::string, Term> terms = compositeResult.getTerms();

	CommutingSetGenerator gen;
	auto commutingSets = gen.getCommutingSet(compositeResult, nQubits);
	auto pi = boost::math::constants::pi<double>();
	auto gateRegistry = xacc::getService<IRProvider>("gate");

	auto uccsdGateFunction = gateRegistry->createFunction("uccsdPrep",{},
			variables);


	// Perform Trotterization...
	for (auto s : commutingSets) {

		for (auto inst : s) {
			Term spinInst = inst;

			// Get the individual pauli terms
			auto termsMap = std::get<2>(spinInst);

			std::vector<std::pair<int,std::string>> terms;
			for (auto& kv : termsMap) {
				if (kv.second != "I" && !kv.second.empty()) {
					terms.push_back({kv.first, kv.second});
				}
			}
			// The largest qubit index is on the last term
			int largestQbitIdx = terms[terms.size() - 1].first;
			auto tempFunction = gateRegistry->createFunction("temp", {}, {});

			for (int i = 0; i < terms.size(); i++) {

				auto qbitIdx = terms[i].first;
				auto gateName = terms[i].second;

				if (i < terms.size() - 1) {
					auto cnot =
							gateRegistry->createInstruction(
									"CNOT",
									std::vector<int> { qbitIdx,
											terms[i + 1].first });
					tempFunction->addInstruction(cnot);
				}

				if (gateName == "X") {
					auto hadamard =
							gateRegistry->createInstruction(
									"H", std::vector<int> { qbitIdx });
					tempFunction->insertInstruction(0, hadamard);
				} else if (gateName == "Y") {
					auto rx =
							gateRegistry->createInstruction(
									"Rx", std::vector<int> { qbitIdx });
					InstructionParameter p(pi / 2.0);
					rx->setParameter(0, p);
					tempFunction->insertInstruction(0, rx);
				}

				// Add the Rotation for the last term
				if (i == terms.size() - 1) {
					// FIXME DONT FORGET DIVIDE BY 2
					std::stringstream ss;
					ss << 2*std::imag(std::get<0>(spinInst)) << " * "
							<< std::get<1>(spinInst);
					auto rz =
							gateRegistry->createInstruction(
									"Rz", std::vector<int> { qbitIdx });

					InstructionParameter p(ss.str());
					rz->setParameter(0, p);
					tempFunction->addInstruction(rz);

				}
			}

			int counter = tempFunction->nInstructions();
			// Add the instruction on the backend of the circuit
			for (int i = terms.size() - 1; i >= 0; i--) {

				auto qbitIdx = terms[i].first;
				auto gateName = terms[i].second;

				if (i < terms.size() - 1) {
					auto cnot =
							gateRegistry->createInstruction(
									"CNOT",
									std::vector<int> { qbitIdx,
											terms[i + 1].first });
					tempFunction->insertInstruction(counter, cnot);
					counter++;
				}

				if (gateName == "X") {
					auto hadamard =
							gateRegistry->createInstruction(
									"H", std::vector<int> { qbitIdx });
					tempFunction->addInstruction(hadamard);
				} else if (gateName == "Y") {
					auto rx =
							gateRegistry->createInstruction(
									"Rx", std::vector<int> { qbitIdx });
					InstructionParameter p(4 * pi - (pi /2.0));
					rx->setParameter(0, p);
					tempFunction->addInstruction(rx);
				}

			}
			// Add to the total UCCSD State Prep function
			for (auto inst : tempFunction->getInstructions()) {
				uccsdGateFunction->addInstruction(inst);
			}
		}
	}

	for (int i = nElectrons-1; i >= 0; i--) {
		auto xGate = gateRegistry->createInstruction(
				"X", std::vector<int>{i});
		uccsdGateFunction->insertInstruction(0,xGate);
	}

	return uccsdGateFunction;
}

}
}

