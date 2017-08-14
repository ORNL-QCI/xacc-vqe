#include "JordanWignerIRTransformation.hpp"
#include "GateFunction.hpp"
#include <boost/math/constants/constants.hpp>

namespace xacc {
namespace vqe {

std::shared_ptr<IR> JordanWignerIRTransformation::transform(
		std::shared_ptr<IR> ir) {

	// We assume we have a FermionIR instance, which contains
	// one FermionKernel, which contains N FermionInstructions, one
	// for each term in the hamiltonian.

	// We want to map that to a Hamiltonian composed of pauli matrices
	// But, we want each term of that to be a separate IR Function.

	// Create a new GateQIR to hold the spin based terms
	auto newIr = std::make_shared<xacc::quantum::GateQIR>();

	auto fermiKernel = ir->getKernels()[0];

	int counter = 0;
	CompositeSpinInstruction total;

	result.clear();

	// Loop over all Fermionic terms...
	for (auto f : fermiKernel->getInstructions()) {

		auto fermionInst = std::dynamic_pointer_cast<FermionInstruction>(f);

		// Get the creation or annihilation sites
		auto termSites = fermionInst->bits();

		// Get the params indicating if termSite is creation or annihilation
		auto params = fermionInst->getParameters();

		auto fermionCoeff = fermionInst->coefficient;
		auto fermionVar = fermionInst->variable;

		CompositeSpinInstruction current;
		for (int i = 0; i < termSites.size(); i++) {
			// Create Zi's
			std::vector<std::pair<int, std::string>> zs;
			for (int j = 0; j <= termSites[i]-1; j++) {
				zs.push_back({j, "Z"});
			}

			SpinInstruction zSpins(zs);
			CompositeSpinInstruction sigPlusMinus;
			auto xi = std::make_shared<SpinInstruction>(
					std::vector<std::pair<int, std::string>> { { termSites[i],
							"X" } });
			sigPlusMinus.addInstruction(xi);

			if (boost::get<int>(params[i])) { // If this is a creation operator
				// Create sigmaPlus
				auto iyi = std::make_shared<SpinInstruction>(
						std::vector<std::pair<int, std::string>> { {
								termSites[i], "Y" } },
						std::complex<double>(0, -1));
				sigPlusMinus.addInstruction(iyi);
			} else {
				auto negiyi = std::make_shared<SpinInstruction>(
						std::vector<std::pair<int, std::string>> { {
								termSites[i], "Y" } },
						std::complex<double>(0, 1));
				sigPlusMinus.addInstruction(negiyi);
			}

			auto temp = 0.5 * sigPlusMinus;
			if (i == termSites.size() - 1) {
				zSpins.variable = fermionVar;
			}
			temp = zSpins * temp;

			if (i == 0) {
				current = temp;
			} else {
				current = current * temp;
			}
		}

		// Case - we could have a current that has no
		// instructions in it. If so add an Identity to it
		if (current.nInstructions() == 0) {
			auto nullInst = std::make_shared<SpinInstruction>(
					std::vector<std::pair<int, std::string>> { { 0, "I" } });
			current.addInstruction(nullInst);
		}

		auto temp = fermionCoeff * current;
		result = result + temp;

	}

	auto resultsStr = result.toString("");
	boost::replace_all(resultsStr, "+", "+\n");
	std::cout << "Jordan Wigner Transformed Fermion to Spin:\nBEGIN\n" << resultsStr << "\nEND\n\n";
	auto pi = boost::math::constants::pi<double>();

	// Populate GateQIR now...
	for (auto inst : result.getInstructions()) {

		// Cast to a Spin Instruction
		auto spinInst = std::dynamic_pointer_cast<SpinInstruction>(inst);

		if (std::fabs(std::real(spinInst->coefficient)) > 1e-9) {
//			std::cout << "Adding " << spinInst->toString("") << "\n";
			// Create a GateFunction and specify that it has
			// a parameter that is the Spin Instruction coefficient
			// that will help us get it to the user for their purposes.
			auto gateFunction = std::make_shared<xacc::quantum::GateFunction>(
					"term" + std::to_string(counter),
					std::vector<InstructionParameter> { InstructionParameter(
							spinInst->coefficient) });

			// Loop over all terms in the Spin Instruction
			// and create instructions to run on the Gate QPU.
			std::vector<std::shared_ptr<xacc::quantum::GateInstruction>> measurements;
			auto terms = spinInst->getTerms();
			for (int i = terms.size() - 1; i >= 0; i--) {
				auto qbit = terms[i].first;
				auto gateName = terms[i].second;
				auto gateRegistry =
						xacc::quantum::GateInstructionRegistry::instance();
				auto meas = gateRegistry->create("Measure", std::vector<int> {
						qbit });
				xacc::InstructionParameter classicalIdx(qbit);
				meas->setParameter(0, classicalIdx);
				measurements.push_back(meas);

				if (gateName == "X") {
					auto hadamard = gateRegistry->create("H", std::vector<int> {
							qbit });
					gateFunction->addInstruction(hadamard);
				} else if (gateName == "Y") {
					auto rx = gateRegistry->create("Rx",
							std::vector<int> { qbit });
					InstructionParameter p(pi / -2.0);
					rx->setParameter(0, p);
					gateFunction->addInstruction(rx);
				}

			}

			for (auto m : measurements) {
				gateFunction->addInstruction(m);
			}
			newIr->addKernel(gateFunction);
			counter++;
		}
	}

	return newIr;
}

CompositeSpinInstruction JordanWignerIRTransformation::getResult() {
	return result;
}

}
}

