#include "JordanWignerIRTransformation.hpp"
#include "XACC.hpp"

namespace xacc {
namespace vqe {

PauliOperator JordanWignerIRTransformation::transform(FermionKernel& kernel) {
	int myStart = 0;
	int myEnd = kernel.nInstructions();

	result.clear();

	fermionKernel = std::make_shared<FermionKernel>(kernel);

	auto instructions = kernel.getInstructions();
	auto instVec = std::vector<InstPtr>(instructions.begin(), instructions.end());

	auto start = std::clock();

	// Loop over all Fermionic terms...
	for (int z = myStart; z < myEnd; ++z) {

		auto f = instVec[z];

		// Get the creation or annihilation sites
		auto termSites = f->bits();

		// Get the params indicating if termSite is creation or annihilation
		auto params = f->getParameters();

		auto coeff = boost::get<std::complex<double>>(params[f->nParameters() - 2]);
		auto fermionVar = boost::get<std::string>(params[f->nParameters() - 1]);

		PauliOperator current(coeff, fermionVar);
		for (int i = 0; i < termSites.size(); i++) {
			std::map<int, std::string> zs;
			auto isCreation = boost::get<int>(params[i]);

			int index = termSites[i];

			std::complex<double> ycoeff =
					isCreation ?
							std::complex<double>(0, -.5) :
							std::complex<double>(0, .5), xcoeff(.5, 0);

			for (int j = 0; j < index; j++) zs.emplace(std::make_pair(j,"Z"));

			current *= PauliOperator(zs)
					* (PauliOperator( { { index, "X" } }, xcoeff)
							+ PauliOperator( { { index, "Y" } }, ycoeff));
		}

		result += current;
	}

//	std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << "\n";

	return result;
}

std::shared_ptr<IR> JordanWignerIRTransformation::transform(
		std::shared_ptr<IR> ir) {
	auto fermiKernel = ir->getKernels()[0];
	return transform(*std::dynamic_pointer_cast<FermionKernel>(fermiKernel)).toXACCIR();
}

}
}

