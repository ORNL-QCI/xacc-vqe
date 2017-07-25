#include "JordanWignerIRTransformation.hpp"
#include "GateQIR.hpp"

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

	int counter = 1;
	for (auto fermionInst : fermiKernel->getInstructions()) {

		auto gateFunction = std::make_shared<xacc::quantum::GateFunction>(
				"term" + std::to_string(counter));

		auto termSites = fermionInst->bits();
		auto params = fermionInst->getParameters();

		for (int i = 0; i < termSites.size(); i++) {
			std::cout << "ADDING SYMBOL " << termSites[i] << ", " << params[i] << "\n";
			// here symbol 'd' is a creation operator
			// and symbol 'a' is an anhililation operator
			if (boost::get<int>(params[i])) { // If this is a creation operator

			} else {

			}

		}


	}

	return newIr;
}

}
}

