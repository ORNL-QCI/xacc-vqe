#include "BravyiKitaevIRTransformation.hpp"
#include "XACC.hpp"
#include "Fenwick.hpp"

namespace xacc {
namespace vqe {


struct CompositeSpinInstruction addBKResults(struct CompositeSpinInstruction x, struct CompositeSpinInstruction y) {
  return x+y;
}

#pragma omp declare reduction( + : CompositeSpinInstruction : \
        				omp_out = addBKResults(omp_out, omp_in)) \
					initializer (omp_priv(omp_orig))

std::shared_ptr<IR> BravyiKitaevIRTransformation::transform(
		std::shared_ptr<IR> ir) {

	// We assume we have a FermionIR instance, which contains
	// one FermionKernel, which contains N FermionInstructions, one
	// for each term in the hamiltonian.

	// We want to map that to a Hamiltonian composed of pauli matrices
	// But, we want each term of that to be a separate IR Function.

	auto fermiKernel = ir->getKernels()[0];

	CompositeSpinInstruction total;
	result.clear();

	int nQubits = std::stoi(xacc::getOption("n-qubits"));

	FenwickTree tree(nQubits);

	// Loop over all Fermionic terms...
	for (int z = 0; z < fermiKernel->nInstructions(); ++z) {
		auto f = fermiKernel->getInstruction(z);

		auto fermionInst = std::dynamic_pointer_cast<FermionInstruction>(f);

		// Get the creation or annihilation sites
		auto termSites = fermionInst->bits();

		// Get the params indicating if termSite is creation or annihilation
		auto params = fermionInst->getParameters();

		auto fermionCoeff = std::complex<double>(std::real(fermionInst->coefficient), 0.0);
		auto fermionVar = fermionInst->variable;

		CompositeSpinInstruction ladderProduct;
		for (int i = 0; i < termSites.size(); i++) {

			auto creationOrAnnihilation = boost::get<int>(params[i]);

			auto index = termSites[i];

			auto paritySet = tree.getParitySet(index);
			auto ancestors = tree.getUpdateSet(index);
			auto ancestorChildren = tree.getRemainderSet(index);

			std::complex<double> dcoeff;
			if (creationOrAnnihilation) {
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

			SpinInstruction d_majorana(dTerms, dcoeff), c_majorana(cTerms, std::complex<double>(0.5,0.0));
			if (i == termSites.size() - 1) {
				d_majorana.variable = fermionVar;
				c_majorana.variable = fermionVar;
			}
			auto sum = c_majorana + d_majorana;
			if (i == 0) {
				ladderProduct = sum;
			} else {
				ladderProduct = ladderProduct * sum;
			}
		}

		// Case - we could have a current that has no
		// instructions in it. If so add an Identity to it
		if (ladderProduct.nInstructions() == 0) {
			auto nullInst = std::make_shared<SpinInstruction>(
					std::vector<std::pair<int, std::string>> { { 0, "I" } });
			ladderProduct.addInstruction(nullInst);
		}

		auto tmp = fermionCoeff * ladderProduct;

		total = total + tmp;
	}

	total.simplify();

	result = total;

	return generateIR();
}

}
}

