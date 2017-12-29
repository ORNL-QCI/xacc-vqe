#include "BravyiKitaevIRTransformation.hpp"
#include "XACC.hpp"
#include "Fenwick.hpp"

namespace xacc {
namespace vqe {

struct CompositeSpinInstruction addBKResults(struct CompositeSpinInstruction x,
		struct CompositeSpinInstruction y) {
	return x + y;
}

#pragma omp declare reduction( + : CompositeSpinInstruction : \
        				omp_out = addBKResults(omp_out, omp_in)) \
					initializer (omp_priv(omp_orig))

std::shared_ptr<IR> BravyiKitaevIRTransformation::transform(
		std::shared_ptr<IR> ir) {

	boost::mpi::communicator world;

	auto fermiKernel = ir->getKernels()[0];

	CompositeSpinInstruction total;
	result.clear();

	int nQubits = std::stoi(xacc::getOption("n-qubits"));

	FenwickTree tree(nQubits);

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

		CompositeSpinInstruction ladderProduct;
		auto nullInst = std::make_shared<SpinInstruction>(
				std::map<int, std::string> { },
				fermionInst->coefficient, fermionVar);
		ladderProduct.addInstruction(nullInst);

		for (int i = 0; i < termSites.size(); i++) {

			auto isCreation = boost::get<int>(params[i]);

			auto index = termSites[i];

			auto paritySet = tree.getParitySet(index);
			auto ancestors = tree.getUpdateSet(index);
			auto ancestorChildren = tree.getRemainderSet(index);

			std::complex<double> dcoeff, ccoeff(.5,0);
			if (isCreation) {
				dcoeff = std::complex<double>(0,-.5);
			} else {
				dcoeff = std::complex<double>(0,.5);
			}

			std::vector<std::pair<int, std::string>> dTerms{{index, "Y"}}, cTerms{{index, "X"}};
			for (auto ac : ancestorChildren) {
				dTerms.push_back({ac->index, "Z"});
			}
			for (auto p : paritySet) {
				cTerms.push_back({p->index, "Z"});
			}
			for (auto a : ancestors) {
				dTerms.push_back({a->index, "X"});
				cTerms.push_back({a->index, "X"});
			}

			SpinInstruction d({},dcoeff), c({},ccoeff);
			for (auto t : dTerms) {
				d = d * SpinInstruction({{t.first,t.second}});
			}
			for (auto t : cTerms) {
				c = c * SpinInstruction({{t.first,t.second}});
			}

			auto sum = c + d;
			ladderProduct = ladderProduct * sum;
			ladderProduct.simplify();
		}

		total = total + ladderProduct;
		total.simplify();
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

