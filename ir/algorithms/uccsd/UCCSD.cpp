#include "GateFunction.hpp"
#include "RuntimeOptions.hpp"
#include "FermionToSpinTransformation.hpp"
#include "ServiceRegistry.hpp"
#include "CommutingSetGenerator.hpp"
#include <boost/math/constants/constants.hpp>
#include "UCCSD.hpp"

using namespace xacc::quantum;

namespace xacc {
namespace vqe {


std::shared_ptr<Function> UCCSD::generate(
		std::shared_ptr<AcceleratorBuffer> buffer,
		std::vector<InstructionParameter> parameters) {

	auto runtimeOptions = RuntimeOptions::instance();

	if (!runtimeOptions->exists("n-electrons")) {
		XACCError("To use this UCCSD State Prep IRGenerator, you "
				"must specify the number of electrons.");
	}

	if (!runtimeOptions->exists("n-qubits")) {
		XACCError("To use this UCCSD State Prep IRGenerator, you "
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
	XACCInfo("Constructing UCCSD Fermion Operator.");
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
				fermiInstruction2->coefficient = -1.0
						* fermiInstruction2->coefficient;
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

							fermiInstruction2->coefficient = -1.0
									* fermiInstruction2->coefficient;
							kernel->addInstruction(fermiInstruction2);

						}
					}
				}
			}
		}
	}

	//std::cout << "KERNEL: \n" << kernel->toString("") << "\n";
	// Create the FermionIR to pass to our transformation.
	auto fermionir = std::make_shared<FermionIR>();
	fermionir->addKernel(kernel);

	XACCInfo("Done constructing UCCSD Fermion Operator.");
	XACCInfo("Mapping UCCSD Fermion Operator to Spin. ");

	std::shared_ptr<IRTransformation> transform;
	if (runtimeOptions->exists("fermion-transformation")) {
		auto transformStr = (*runtimeOptions)["fermion-transformation"];
		transform = ServiceRegistry::instance()->getService<IRTransformation>(
				transformStr);
	} else {
		transform = ServiceRegistry::instance()->getService<IRTransformation>(
				"jordan-wigner");
	}

	XACCInfo("Done mapping UCCSD Fermion Operator to Spin. ");

	std::dynamic_pointer_cast<FermionToSpinTransformation>(transform)->runParallel = false;

	// Create the Spin Hamiltonian
	auto transformedIR = transform->transform(fermionir);

	auto compositeResult =
			std::dynamic_pointer_cast<FermionToSpinTransformation>(transform)->getResult();

	// Convert imaginary part to real part
	for (auto term : compositeResult.getInstructions()) {
		auto castedTerm = std::dynamic_pointer_cast<SpinInstruction>(term);
		castedTerm->coefficient = std::complex<double>(
				std::imag(castedTerm->coefficient), 0);
	}

	CommutingSetGenerator gen;
	auto commutingSets = gen.getCommutingSet(compositeResult, nQubits);
	auto pi = boost::math::constants::pi<double>();
	auto uccsdGateFunction = std::make_shared<xacc::quantum::GateFunction>("uccsdPrep",
			variables);

	// Perform Trotterization...
	for (auto s : commutingSets) {

		for (auto instIdx : s) {
			auto temp = compositeResult.getInstruction(instIdx);
			auto spinInst = std::dynamic_pointer_cast<SpinInstruction>(temp);
			// Get the individual pauli terms
			auto terms = spinInst->getTerms();

			// The largest qubit index is on the last term
			int largestQbitIdx = terms[terms.size() - 1].first;
			auto tempFunction = std::make_shared<xacc::quantum::GateFunction>(
					"");

			for (int i = 0; i < terms.size(); i++) {

				auto qbitIdx = terms[i].first;
				auto gateName = terms[i].second;

				if (i < terms.size() - 1) {
					auto cnot =
							xacc::quantum::GateInstructionRegistry::instance()->create(
									"CNOT",
									std::vector<int> { qbitIdx,
											terms[i + 1].first });
					tempFunction->addInstruction(cnot);
				}

				if (gateName == "X") {
					auto hadamard =
							xacc::quantum::GateInstructionRegistry::instance()->create(
									"H", std::vector<int> { qbitIdx });
					tempFunction->insertInstruction(0, hadamard);
				} else if (gateName == "Y") {
					auto rx =
							xacc::quantum::GateInstructionRegistry::instance()->create(
									"Rx", std::vector<int> { qbitIdx });
					InstructionParameter p(pi / 2.0);
					rx->setParameter(0, p);
					tempFunction->insertInstruction(0, rx);
				}

				// Add the Rotation for the last term
				if (i == terms.size() - 1) {
					// FIXME DONT FORGET DIVIDE BY 2
					std::stringstream ss;
					ss << 2*std::real(spinInst->coefficient) << " * "
							<< spinInst->variable;
					auto rz =
							xacc::quantum::GateInstructionRegistry::instance()->create(
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
							xacc::quantum::GateInstructionRegistry::instance()->create(
									"CNOT",
									std::vector<int> { qbitIdx,
											terms[i + 1].first });
					tempFunction->insertInstruction(counter, cnot);
					counter++;
				}

				if (gateName == "X") {
					auto hadamard =
							xacc::quantum::GateInstructionRegistry::instance()->create(
									"H", std::vector<int> { qbitIdx });
					tempFunction->addInstruction(hadamard);
				} else if (gateName == "Y") {
					auto rx =
							xacc::quantum::GateInstructionRegistry::instance()->create(
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
		auto xGate = xacc::quantum::GateInstructionRegistry::instance()->create(
				"X", std::vector<int>{i});
		uccsdGateFunction->insertInstruction(0,xGate);
	}

	return uccsdGateFunction;
}

}
}

