#include "JordanWignerIRTransformation.hpp"
#include "GateFunction.hpp"

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

		auto temp = fermionCoeff * current;
		result = result + temp;

	}

	std::cout << "Transformed: " << result.toString("") << "\n";
	// Populate GateQIR now...
	for (auto inst : result.getInstructions()) {

		// Cast to a Spin Instruction
		auto spinInst = std::dynamic_pointer_cast<SpinInstruction>(inst);

		// Create a GateFunction and specify that it has
		// a parameter that is the Spin Instruction coefficient
		// that will help us get it to the user for their purposes.
		auto gateFunction = std::make_shared<xacc::quantum::GateFunction>(
				"term" + std::to_string(counter),
				std::vector<InstructionParameter> { InstructionParameter(
						spinInst->coefficient) });

		// Loop over all terms in the Spin Instruction
		// and create instructions to run on the Gate QPU.
		auto terms = spinInst->getTerms();
		for (int i = terms.size()-1; i >= 0; i--) {
			auto qbit = terms[i].first;
			auto gateName = terms[i].second;
			std::shared_ptr<xacc::quantum::GateInstruction> meas;
			auto gateRegistry = xacc::quantum::GateInstructionRegistry::instance();

			if (gateName != "I") {
				meas = gateRegistry->create("Measure", std::vector<int>{qbit});
				xacc::InstructionParameter classicalIdx(qbit);
				meas->setParameter(0, classicalIdx);
			} else {
				// FIXME IF IT IS I, FIGURE THIS OUT....
			}

			// If its an X we have to add a Hadamard before the measure
			// If its a Y we have to add a Rx(-pi / 2) gate before the measure.
			if (gateName == "X") {
				gateFunction->addInstruction(gateRegistry->create("H", std::vector<int>{qbit}));
			} else if (gateName == "Y") {
				auto rx = gateRegistry->create("Rx", std::vector<int>{qbit});
				InstructionParameter p(3.1415926 / -2.0);
				rx->setParameter(0, p);
				gateFunction->addInstruction(rx);
			}

			// Add the MeasureZ gate...
			gateFunction->addInstruction(meas);
		}

		newIr->addKernel(gateFunction);
		counter++;
	}

	return newIr;
}

CompositeSpinInstruction JordanWignerIRTransformation::getResult() {
	return result;
}

}
}

