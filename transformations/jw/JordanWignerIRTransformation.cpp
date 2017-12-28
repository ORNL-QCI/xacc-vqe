#include "JordanWignerIRTransformation.hpp"
#include "XACC.hpp"


namespace xacc {
namespace vqe {

struct CompositeSpinInstruction addJWResults(struct CompositeSpinInstruction x,
		struct CompositeSpinInstruction y) {
	return x + y;
}

#pragma omp declare reduction( + : CompositeSpinInstruction : \
		omp_out = addJWResults(omp_out, omp_in)) \
		initializer( omp_priv(omp_orig))

std::shared_ptr<IR> JordanWignerIRTransformation::transform(
		std::shared_ptr<IR> ir) {

	boost::mpi::communicator world;

	auto fermiKernel = ir->getKernels()[0];

	CompositeSpinInstruction total;
	result.clear();

	int myStart = 0;
	int myEnd = fermiKernel->nInstructions();

	if (runParallel) {
		myStart = (world.rank()) * fermiKernel->nInstructions() / world.size();
		myEnd = (world.rank() + 1) * fermiKernel->nInstructions() / world.size();
	}

	// Loop over all Fermionic terms...
#pragma omp parallel for shared(fermiKernel) reduction (+:total) if (runParallel)
	for (int z = myStart; z < myEnd; ++z) {

		auto f = fermiKernel->getInstruction(z);

		auto fermionInst = std::dynamic_pointer_cast<FermionInstruction>(f);

		// Get the creation or annihilation sites
		auto termSites = fermionInst->bits();

		// Get the params indicating if termSite is creation or annihilation
		auto params = fermionInst->getParameters();

		auto fermionVar = fermionInst->variable;

		CompositeSpinInstruction current;
		auto nullInst = std::make_shared<SpinInstruction>(
				std::vector<std::pair<int, std::string>> { { 0, "I" } },
				fermionInst->coefficient, fermionVar);
		current.addInstruction(nullInst);

		for (int i = 0; i < termSites.size(); i++) {
			auto isCreation = boost::get<int>(params[i]);

			int index = termSites[i];

			std::complex<double> ycoeff, xcoeff(.5,0);
			if (isCreation) {
				ycoeff = std::complex<double>(0,-.5);
			} else {
				ycoeff = std::complex<double>(0,.5);
			}

			// Create Zi's
			std::vector<std::pair<int, std::string>> zx, zy;
			for (int j = 0; j < index; j++) {
				zx.push_back({j, "Z"});
				zy.push_back({j, "Z"});
			}

			zx.push_back({index, "X"});
			zy.push_back({index, "Y"});

			SpinInstruction xcomp(zx, xcoeff), ycomp(zy, ycoeff);
			auto sum = xcomp + ycomp;
			current = current * sum;
			current.simplify();
//			current.compress();
		}

		total = total + current;
		total.simplify();
//		total.compress();
	}

	total.simplify();

	CompositeSpinInstruction i;

	if (world.size() > 1 && runParallel) {
		i = distributeResults(world, total);
	} else {
		i = total;
	}

	world.barrier();

	result = i;

	result.compress();

	return generateIR();
}

}
}

