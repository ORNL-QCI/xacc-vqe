#include "GateFunction.hpp"
#include <boost/math/constants/constants.hpp>
#include "BravyiKitaevIRTransformation.hpp"
#include "XACC.hpp"
#include "Fenwick.hpp"

namespace xacc {
namespace vqe {

std::shared_ptr<IR> BravyiKitaevIRTransformation::transform(
		std::shared_ptr<IR> ir) {

	// We assume we have a FermionIR instance, which contains
	// one FermionKernel, which contains N FermionInstructions, one
	// for each term in the hamiltonian.

	// We want to map that to a Hamiltonian composed of pauli matrices
	// But, we want each term of that to be a separate IR Function.

	// Create a new GateQIR to hold the spin based terms
	auto newIr = std::make_shared<xacc::quantum::GateQIR>();

	auto fermiKernel = ir->getKernels()[0];

	int counter = 0;
	CompositeSpinInstruction total;
	result.clear();

	int nQubits = std::stoi(xacc::getOption("n-qubits"));

	FenwickTree tree(nQubits);


	// Loop over all Fermionic terms...
	for (auto f : fermiKernel->getInstructions()) {

		auto fermionInst = std::dynamic_pointer_cast<FermionInstruction>(f);

		// Get the creation or annihilation sites
		auto termSites = fermionInst->bits();

		// Get the params indicating if termSite is creation or annihilation
		auto params = fermionInst->getParameters();

		auto fermionCoeff = fermionInst->coefficient;
		auto fermionVar = fermionInst->variable;

		CompositeSpinInstruction ladderProduct;
		for (int i = 0; i < termSites.size(); i++) {

			auto creationOrAnnihilation = boost::get<int>(params[i]);

			auto index = termSites[i];

			auto paritySet = tree.getParitySet(index);
			auto ancestors = tree.getUpdateSet(index);
			auto ancestorChildren = tree.getRemainderSet(index);

			auto getSetElement = [](const std::set<std::shared_ptr<Node>>& s, const int idx) {
				return *std::next(s.begin(), idx);
			};

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

		ladderProduct = fermionCoeff * ladderProduct;

		total = total + ladderProduct;

	}

	result = total;

	auto resultsStr = result.toString("");
	boost::replace_all(resultsStr, "+", "+\n");
	std::cout << "BK Transformed Fermion to Spin:\nBEGIN\n" << resultsStr << "\nEND\n\n";
	auto pi = boost::math::constants::pi<double>();

	// Populate GateQIR now...
	for (auto inst : result.getInstructions()) {

		// Cast to a Spin Instruction
		auto spinInst = std::dynamic_pointer_cast<SpinInstruction>(inst);

		if (std::fabs(std::real(spinInst->coefficient)) > 1e-9) {

			int isIdentity = 0;
			if (*spinInst.get() == SpinInstruction( { { 0, "I" } })) {
				isIdentity = 1;
			}

//			std::cout << "Adding " << spinInst->toString("") << "\n";
			// Create a GateFunction and specify that it has
			// a parameter that is the Spin Instruction coefficient
			// that will help us get it to the user for their purposes.
			auto gateFunction = std::make_shared<xacc::quantum::GateFunction>(
					"term" + std::to_string(counter),
					std::vector<InstructionParameter> { InstructionParameter(
							spinInst->coefficient), InstructionParameter(
									isIdentity) });

			// Loop over all terms in the Spin Instruction
			// and create instructions to run on the Gate QPU.
			std::vector<std::shared_ptr<xacc::quantum::GateInstruction>> measurements;
			auto terms = spinInst->getTerms();
			for (int i = terms.size() - 1; i >= 0; i--) {
				auto qbit = terms[i].first;
				auto gateName = terms[i].second;
				auto gateRegistry =
						xacc::quantum::GateInstructionRegistry::instance();
				auto meas = gateRegistry->create("Measure", std::vector<int> {
						qbit });
				xacc::InstructionParameter classicalIdx(qbit);
				meas->setParameter(0, classicalIdx);
				measurements.push_back(meas);

				if (gateName == "X") {
					auto hadamard = gateRegistry->create("H", std::vector<int> {
							qbit });
					gateFunction->addInstruction(hadamard);
				} else if (gateName == "Y") {
					auto rx = gateRegistry->create("Rx",
							std::vector<int> { qbit });
					InstructionParameter p(pi / 2.0);
					rx->setParameter(0, p);
					gateFunction->addInstruction(rx);
				}

			}

			for (auto m : measurements) {
				gateFunction->addInstruction(m);
			}
			newIr->addKernel(gateFunction);
			counter++;
		}
	}


	return newIr;
}

CompositeSpinInstruction BravyiKitaevIRTransformation::getResult() {
	return result;
}

}
}

