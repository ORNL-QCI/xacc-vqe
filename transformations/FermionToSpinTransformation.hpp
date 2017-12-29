/***********************************************************************************
 * Copyright (c) 2017, UT-Battelle
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of the xacc nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Contributors:
 *   Initial API and implementation - Alex McCaskey
 *
 **********************************************************************************/
#ifndef VQE_IR_FERMIONTOSPINTRANSFORMATION_HPP_
#define VQE_IR_FERMIONTOSPINTRANSFORMATION_HPP_

#include "IRTransformation.hpp"
#include "FermionKernel.hpp"
#include "FermionIR.hpp"
#include "GateQIR.hpp"
#include "SpinInstruction.hpp"
#include "GateFunction.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/mpi.hpp>

#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include "XACC.hpp"

namespace xacc {

namespace vqe {

using PersistedSpinInstruction = std::pair<std::complex<double>, std::vector<int>>;

struct add_spin_instructions {
	std::vector<PersistedSpinInstruction> operator()( std::vector<PersistedSpinInstruction> a, std::vector<PersistedSpinInstruction> b) {
		std::move(b.begin(), b.end(), std::back_inserter(a));
		int nQubits = std::stoi(xacc::getOption("n-qubits"));
		CompositeSpinInstruction i;
		std::vector<PersistedSpinInstruction> persistedInstructions;
		for (auto inst : a) {
			auto spin = std::make_shared<SpinInstruction>();
			spin->fromBinaryVector(inst.second, inst.first);
			i.addInstruction(spin);
		}

		i.simplify();

		for (auto inst : i.getInstructions()) {
			auto casted = std::dynamic_pointer_cast<SpinInstruction>(inst);
			auto x = casted->toBinaryVector(nQubits);
			persistedInstructions.push_back({casted->coefficient, x});
		}

		return persistedInstructions;
	}
};
/**
 */
class FermionToSpinTransformation: public xacc::IRTransformation {

public:

	/**
	 * Return the result of the Jordan Wigner transformation
	 * @return
	 */
	virtual CompositeSpinInstruction getResult() {
		return result;
	}

	bool runParallel = true;

protected:

	/**
	 * Reference to the transformation result.
	 */
	CompositeSpinInstruction result;

	CompositeSpinInstruction distributeResults(boost::mpi::communicator& world, CompositeSpinInstruction total) {
		CompositeSpinInstruction i;
		std::vector<PersistedSpinInstruction> persistedInstructions, globalInstructions;
		int nQubits = std::stoi(xacc::getOption("n-qubits"));
		for (auto inst : total.getInstructions()) {
			auto casted = std::dynamic_pointer_cast<SpinInstruction>(inst);
			std::vector<int> x = casted->toBinaryVector(nQubits);
			persistedInstructions.push_back({casted->coefficient, x});
		}
		boost::mpi::all_reduce(world, persistedInstructions, globalInstructions, add_spin_instructions());
		for (auto inst : globalInstructions) {
			auto spin = std::make_shared<SpinInstruction>();
			spin->fromBinaryVector(inst.second, inst.first);
			i.addInstruction(spin);
		}
		return i;
	}

	std::shared_ptr<IR> generateIR() {
		boost::mpi::communicator world;

		// Create a new GateQIR to hold the spin based terms
		auto newIr = std::make_shared<xacc::quantum::GateQIR>();
		int counter = 0;
		auto resultsStr = result.toString("");
		boost::replace_all(resultsStr, "+", "+\n");
		if (world.rank() == 0) std::cout << "Transformed Fermion to Spin:\nBEGIN\n" << resultsStr << "\nEND\n\n";
		auto pi = boost::math::constants::pi<double>();
//		exit(0);
		// Populate GateQIR now...
		for (auto inst : result.getInstructions()) {

			// Cast to a Spin Instruction
			auto spinInst = std::dynamic_pointer_cast<SpinInstruction>(inst);

//			if (std::fabs(std::real(spinInst->coefficient)) > 1e-9) {

//				int isIdentity = 0;
//				if (spinInst->isIdentity()) {
//					isIdentity = 1;
//				}

				// Create a GateFunction and specify that it has
				// a parameter that is the Spin Instruction coefficient
				// that will help us get it to the user for their purposes.
				auto gateFunction = std::make_shared<xacc::quantum::GateFunction>(
						"term" + std::to_string(counter),
						std::vector<InstructionParameter> { InstructionParameter(
								spinInst->coefficient), InstructionParameter(
										spinInst->isIdentity() ? 1 : 0) });

				// Loop over all terms in the Spin Instruction
				// and create instructions to run on the Gate QPU.
				std::vector<std::shared_ptr<xacc::quantum::GateInstruction>> measurements;
				auto termsMap = spinInst->getTerms();

				std::vector<std::pair<int,std::string>> terms;
				for (auto& kv : termsMap) {
					terms.push_back({kv.first, kv.second});
				}

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

				if (!spinInst->isIdentity()) {
					for (auto m : measurements) {
						gateFunction->addInstruction(m);
					}
				}
				newIr->addKernel(gateFunction);
				counter++;
			}
//		}
		return newIr;
	}
};

}

}
#endif
