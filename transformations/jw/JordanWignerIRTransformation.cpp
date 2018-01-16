#include "JordanWignerIRTransformation.hpp"
#include "XACC.hpp"

namespace xacc {
namespace vqe {

PauliOperator JordanWignerIRTransformation::transform(FermionKernel& kernel) {
	PauliOperator localResult;
	int myStart = 0;
	int myEnd = kernel.nInstructions();

	result.clear();

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
			auto isCreation = boost::get<int>(params[i]);

			int index = termSites[i];

			std::complex<double> ycoeff, xcoeff(.5,0);
			if (isCreation) {
				ycoeff = std::complex<double>(0,-.5);
			} else {
				ycoeff = std::complex<double>(0,.5);
			}

			PauliOperator zs(1.0);
			for (int j = 0; j < index; j++) {
				zs *= PauliOperator({{j, "Z"}});
			}

			current *= zs
					* (PauliOperator( { { index, "X" } }, xcoeff)
							+ PauliOperator( { { index, "Y" } }, ycoeff));
		}

		localResult += current;
	}

	std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << "\n";

	result = localResult;

	return localResult;
}

std::shared_ptr<IR> JordanWignerIRTransformation::transform(
		std::shared_ptr<IR> ir) {
	boost::mpi::communicator world;
	auto fermiKernel = ir->getKernels()[0];
	return transform(*std::dynamic_pointer_cast<FermionKernel>(fermiKernel)).toXACCIR();
}

}
}

