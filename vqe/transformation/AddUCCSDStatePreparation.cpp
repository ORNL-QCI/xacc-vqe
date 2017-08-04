#include "AddUCCSDStatePreparation.hpp"
#include "GateFunction.hpp"
#include "RuntimeOptions.hpp"
#include "FermionToSpinTransformation.hpp"
#include "ServiceRegistry.hpp"

namespace xacc {
namespace vqe {

std::shared_ptr<IR> AddUCCSDStatePreparation::transform(
		std::shared_ptr<IR> ir) {

	auto runtimeOptions = RuntimeOptions::instance();

	// Create a new GateQIR to hold the spin based terms
	auto newIr = std::make_shared<xacc::quantum::GateQIR>();

	if (!runtimeOptions->exists("n-electrons")) {
		XACCError("To use this State Prep Transformation, you "
				"must specify the number of electrons.");
	}

	if (!runtimeOptions->exists("n-qubits")) {
		XACCError("To use this State Prep Transformation, you "
				"must specify the number of qubits.");
	}

	auto nQubits = std::stoi((*runtimeOptions)["n-qubits"]);
	auto nElectrons = std::stoi((*runtimeOptions)["n-electrons"]);

	// Compute the number of parameters
	auto nOccupied = (int) std::ceil(nElectrons / 2.0);
	auto nVirtual = nQubits / 2 - nOccupied;
	auto nSingle = nOccupied * nVirtual;
	auto nDouble = std::pow(nSingle, 2);
	auto nParams = nSingle + nDouble;

	auto singletIndex = [=](int i, int j) -> int {
	        return i * nOccupied + j;
	};

	auto doubletIndex = [=](int i, int j, int k, int l) -> int {
		return
		(i * nOccupied * nVirtual * nOccupied +
				j * nVirtual * nOccupied +
				k * nOccupied +
				l);
	};

	std::vector<std::string> params;
	for (int i = 0; i < nParams; i++) {
		params.push_back("theta"+std::to_string(i));
	}

	std::cout << "HEY: " << nOccupied << ", " << nVirtual << ", " << nParams << "\n";

	auto kernel = std::make_shared<FermionKernel>("fermiUCCSD");

	for (int i = 0; i < nVirtual; i++) {
		for (int j = 0; j < nOccupied; j++) {
			for (int l = 0; l < 2; l++) {
				std::cout << i << " " << j << " " << l << "\n";
				std::vector<std::pair<int, int>> operators { { 2
						* (i + nOccupied) + l, 1 }, { 2 * j + l, 0 } };
				auto fermiInstruction1 = std::make_shared<FermionInstruction>(
						operators, params[singletIndex(i, j)]);
				kernel->addInstruction(fermiInstruction1);

				std::vector<std::pair<int, int>> operators2 { { 2 * j + l, 1}, {2 * (i + nOccupied) + l, 0} };
				auto fermiInstruction2 = std::make_shared<FermionInstruction>(
						operators2, params[singletIndex(i, j)]);
				fermiInstruction2->coefficient = -1.0 * fermiInstruction2->coefficient;
				kernel->addInstruction(fermiInstruction2);

			}
		}
	}

	for (int i = 0; i < nVirtual; i++) {
		for (int j = 0; j < nOccupied; j++) {
			for (int l = 0; l < 2; l++) {
				for (int i2 = 0; i2 < nVirtual; i2++) {
					for (int j2 = 0; j2 < nOccupied; j2++) {
						for (int l2 = 0; l2 < 2; l2++) {
							std::vector<std::pair<int, int>> operators1 { { 2
									* (i + nOccupied) + l, 1 },
									{ 2 * j + l, 0 }, { 2 * (i2 + nOccupied)
											+ l2, 1 }, { 2 * j2 + l2, 0 } };

							std::vector<std::pair<int, int>> operators2 { { 2
									* j2 + l2, 1 }, { 2 * (i2 + nOccupied) + l2,
									0 }, { 2 * j + l, 1 }, { 2 * (i + nOccupied)
									+ l, 0 } };

							auto doubletIdx1 = nSingle + doubletIndex(i, j, i2, j2);
							// FIXME THIS HAS TO BE NEGATIVE
							auto doubletIdx2 = nSingle + doubletIndex(i, j, i2, j2);

							auto fermiInstruction1 = std::make_shared<
									FermionInstruction>(operators1,
									params[doubletIdx1]);

							kernel->addInstruction(fermiInstruction1);

							auto fermiInstruction2 = std::make_shared<
									FermionInstruction>(operators2,
									params[doubletIdx2]);

							fermiInstruction2->coefficient = -1.0 * fermiInstruction2->coefficient;
							kernel->addInstruction(fermiInstruction2);

						}
					}
				}
			}
		}
	}

	// Create the FermionIR to pass to our transformation.
	auto fermionir = std::make_shared<FermionIR>();
	fermionir->addKernel(kernel);

	std::cout << "HEY: " << kernel->toString("") << "\n";

	std::shared_ptr<IRTransformation> transform;
	if (runtimeOptions->exists("fermion-transformation")) {
		auto transformStr = (*runtimeOptions)["fermion-transformation"];
		transform = ServiceRegistry::instance()->getService<IRTransformation>(
				transformStr);
	} else {
		transform = ServiceRegistry::instance()->getService<IRTransformation>(
				"jordan-wigner");
	}

	// Create the Spin Hamiltonian
	auto transformedIR = transform->transform(fermionir);

	auto compositeResult =
			std::dynamic_pointer_cast<FermionToSpinTransformation>(transform)->getResult();

	// Compute commuting sets...
	auto supports = [](CompositeSpinInstruction& i, int idx) -> std::vector<std::pair<int, std::string>> {
		auto spinInst = std::dynamic_pointer_cast<SpinInstruction>(i.getInstruction(idx));
		return spinInst->getTerms();
	};

	auto commutator =
			[](std::vector<std::pair<int, std::string>> support1,
					std::vector<std::pair<int, std::string>> support2) -> bool {
		std::vector<int> iqbits, jqbits;
		std::vector<std::string> ips, jps;

		for (auto t : support1) {
			iqbits.push_back(t.first);
			ips.push_back(t.second);
		}
		for (auto t : support2) {
			jqbits.push_back(t.first);
			jps.push_back(t.second);
		}

		std::vector<int> overlaps;
		for (int i = 0; i < iqbits.size(); i++) {
			auto site = iqbits[i];
			auto itr = std::find(jqbits.begin(), jqbits.end(), site);
			if (itr != jqbits.end()) {
				auto ind2 = std::distance(jqbits.begin(), itr);
				if (ips[i] != jps[ind2]) {
					overlaps.push_back(site);
				}
			}
		}
		return (bool) (overlaps.size() + 1 % 2);
	};
/*
	commuting_sets = [] # collect groups of commuting operators by thier indicies.
	number_of_ops = len(tuple(composite_spin_operator.terms.keys()))

	for i in range(number_of_ops):
	    if i == 0:
	        # first operator initilizes first set
	        commuting_sets.append(set([i]))

	    # check commutators with earlier terms
	    for j in range(i):
	        if commutator(supports(composite_spin_operator, i),
	                      supports(composite_spin_operator, j)) == True:
	            print('found commuting term -- existing set expanded to:')
	            for s in commuting_sets:
	                if j in s:
	                    s.add(i)
	                    print(s)
	            break

	    # if no commutators found
	    if not any([i in cs for cs in commuting_sets]):
	        commuting_sets.append(set([i]))


	commuting_sets */
	std::vector<std::vector<int>> commutingSets;
	for (int i = 0; i < compositeResult.getInstructions().size(); i++) {
		if (i == 0) {
			commutingSets.push_back(std::vector<int>{i});
		}

		for (int j = 0 ; j < i; j++) {
			if (commutator(supports(compositeResult, i), supports(compositeResult,j))) {
				for (auto s : commutingSets) {
					auto itr = std::find(s.begin(), s.end(), j);
					if (itr != s.end()) {
						s.push_back(i);
					}
				}
				break;
			}
		}

		bool found = false;
		for (auto s : commutingSets) {
			auto itr = std::find(s.begin(), s.end(), i);
			if (itr != s.end()) {
				found = true;
				break;
			}
		}
		if (!found) {
			std::cout << "NOT FOUND\n";
			commutingSets.push_back(std::vector<int>{i});
		}
	}

	std::cout << "Commuting term indices:\n";
	for (auto s : commutingSets) {
		std::cout << "{ ";
		for (auto si : s) {
			std::cout << si << " ";
		}
		std::cout << "}";
	}

	// Perform Trotterization...

	// Prepend all GateFunctions with function containing trotterized circuit

	return newIr;
}

}
}

