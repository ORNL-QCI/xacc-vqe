#include "BravyiKitaevIRTransformation.hpp"
#include "XACC.hpp"
#include "Fenwick.hpp"

namespace xacc {
namespace vqe {

PauliOperator BravyiKitaevIRTransformation::transform(FermionKernel& kernel) {
	result.clear();

	int nQubits = std::stoi(xacc::getOption("n-qubits"));

	FenwickTree tree(nQubits);

	int myStart = 0;
	int myEnd = kernel.nInstructions();

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

		PauliOperator ladderProduct(coeff, fermionVar);

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

			std::map<int, std::string> dTerms{{index, "Y"}}, cTerms{{index, "X"}};
			for (auto ac : ancestorChildren) {
				dTerms[ac->index] = "Z";
			}
			for (auto p : paritySet) {
				cTerms[p->index] = "Z";
			}
			for (auto a : ancestors) {
				dTerms[a->index] = "X";
				cTerms[a->index] = "X";
			}

			PauliOperator d(dTerms,dcoeff), c(cTerms,ccoeff);

			ladderProduct *= (c+d);
		}

		result += ladderProduct;
	}
	std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << "\n";

	return result;
}

std::shared_ptr<IR> BravyiKitaevIRTransformation::transform(
		std::shared_ptr<IR> ir) {
	auto fermiKernel = ir->getKernels()[0];
	return transform(*std::dynamic_pointer_cast<FermionKernel>(fermiKernel)).toXACCIR();
}

}
}

