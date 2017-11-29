#include "JordanWignerIRTransformation.hpp"
#include "GateFunction.hpp"
#include <boost/math/constants/constants.hpp>

namespace xacc {
namespace vqe {

struct CompositeSpinInstruction addJWResults(struct CompositeSpinInstruction x, struct CompositeSpinInstruction y) {
  return x+y;
}

#pragma omp declare reduction( + : CompositeSpinInstruction : omp_out = omp_out + omp_in ) \
  initializer( omp_priv = omp_orig )

//#pragma omp declare reduction( CustomAdd : CompositeSpinInstruction : \
//       				omp_out = addJWResults(omp_out, omp_in)) 
//initializer (omp_priv(omp_orig))

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

	CompositeSpinInstruction total;

	result.clear();

	// Loop over all Fermionic terms...
#pragma omp parallel for shared(fermiKernel) reduction (+:total)
	for (int z = 0; z < fermiKernel->nInstructions(); ++z) {

		auto f = fermiKernel->getInstruction(z);

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
		total = total + temp;
	}

	total.simplify();
	result = total;

	return generateIR();
}

}
}

