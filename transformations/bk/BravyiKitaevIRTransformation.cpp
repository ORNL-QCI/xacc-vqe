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
		total.simplify();
	}
	total.simplify();

	CompositeSpinInstruction i;

	if (world.size() > 1 && runParallel) {
		std::string globalString = "";

		// FIXME DO THIS BETTER WITH BINARY VECTOR REPRESENTATION
		auto totalAsString = total.toString("");

		boost::mpi::all_reduce(world, totalAsString, globalString, add_them());

		boost::char_separator<char> sep("+");
		boost::tokenizer<boost::char_separator<char> > tokens(globalString, sep);
		for (auto t : tokens) {
			boost::trim(t);

			std::vector<std::string> splitMult, splitComma;
			boost::split(splitMult, t, boost::is_any_of("*"));

			auto coeffStr = splitMult[0];
			boost::replace_all(coeffStr, "(", "");
			boost::replace_all(coeffStr, ")", "");
			boost::split(splitComma, coeffStr, boost::is_any_of(","));
			auto coeff = std::complex<double>(std::stod(splitComma[0]), std::stod(splitComma[1]));

			std::vector<std::pair<int, std::string>> term;
			for (int i = 1; i < splitMult.size(); i++) {
				auto currentStr = splitMult[i];
				boost::trim(currentStr);

				if (currentStr == "I") {
					term.push_back({0, "I"});
				} else {
					std::stringstream s1, s2;
					s1 << currentStr.at(1);
					s2 << currentStr.at(0);
					auto idxStr = s1.str();
					boost::trim(idxStr);

					int idx = std::stoi(idxStr);
					std::string pauli = s2.str();
					boost::trim(pauli);
					term.push_back( { idx, pauli});
				}
			}

			i.addInstruction(std::make_shared<SpinInstruction>(term, coeff));
		}

		i.simplify();
	} else {
		i = total;
	}

	world.barrier();

	result = i;

	return generateIR();
}

}
}

