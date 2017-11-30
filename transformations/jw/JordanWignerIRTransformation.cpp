#include "JordanWignerIRTransformation.hpp"
#include "GateFunction.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/mpi.hpp>

namespace xacc {
namespace vqe {

using TermType = std::pair<std::complex<double>, std::vector<std::pair<int, std::string>>>;

struct CompositeSpinInstruction addJWResults(struct CompositeSpinInstruction x, struct CompositeSpinInstruction y) {
  return x+y;
}

#pragma omp declare reduction( + : CompositeSpinInstruction : omp_out = omp_out + omp_in ) \
  initializer( omp_priv = omp_orig )

std::shared_ptr<IR> JordanWignerIRTransformation::transform(
		std::shared_ptr<IR> ir) {

	boost::mpi::communicator world;

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

