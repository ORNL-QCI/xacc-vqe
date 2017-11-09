#include "ComputeExpectationValues.hpp"

#include "VQEProgram.hpp"

namespace xacc {
namespace vqe {

VQETaskResult ComputeExpectationValues::execute(
		Eigen::VectorXd parameters) {

	// Local Declarations
	auto comm = program->getCommunicator();
	double sum = 0.0, localExpectationValue = 0.0;
	int rank = comm.rank(), nlocalqpucalls = 0;
	int nRanks = comm.size();

	auto statePrep = program->getStatePreparationCircuit();
	auto nQubits = program->getNQubits();
	auto qpu = program->getAccelerator();

	// Evaluate our variable parameterized State Prep circuite
	// to produce a state prep circuit with actual rotations
	auto evaluatedStatePrep = StatePreparationEvaluator::evaluateCircuit(
			statePrep, program->getNParameters(), parameters);

	auto kernels = program->getVQEKernels();

	VQETaskResult expVals;
	if (qpu->isPhysical()) {

		XACCInfo("Computing Energy with XACC Kernel Multi-Exec.");
		std::vector<double> coeffs;
		KernelList<> modifiedKernelList(qpu);
		KernelList<> identityKernels(qpu);

		for (auto& k : kernels) {
			// If not identity
			if (k.getIRFunction()->nInstructions() > 0) {
				k.getIRFunction()->insertInstruction(0, evaluatedStatePrep);
				modifiedKernelList.push_back(k);
				coeffs.push_back(
						std::real(
								boost::get<std::complex<double>>(
										k.getIRFunction()->getParameter(0))));
			} else {
				identityKernels.push_back(k);
			}
		}

		auto buff = qpu->createBuffer("qreg", nQubits);

		auto tmpBuffers = modifiedKernelList.execute(buff);
		nlocalqpucalls += tmpBuffers.size();

		for (auto k : identityKernels) {
			expVals.push_back({parameters, std::real(
					boost::get<std::complex<double>>(
							k.getIRFunction()->getParameter(0)))});
		}


		int counter = 0;
		for (auto b : tmpBuffers) {
			localExpectationValue = b->getExpectationValueZ();
			expVals.push_back({parameters, coeffs[counter] * localExpectationValue});
			counter++;
		}


		for (int i = 0; i < kernels.size(); i++) {
			if (kernels[i].getIRFunction()->nInstructions() > 0) {
				kernels[i].getIRFunction()->removeInstruction(0);
			}
		}
	} else {
		for (int i = 0; i < kernels.size(); i++) {

			// Get the ith Kernel
			auto kernel = kernels[i];

			if (kernel.getIRFunction()->nInstructions() > 0) {
				// Insert the state preparation circuit IR
				// at location 0 in this Kernels IR instructions.
				kernel.getIRFunction()->insertInstruction(0,
						evaluatedStatePrep);

				// Create a temporary buffer of qubits
				auto buff = qpu->createBuffer("qreg", nQubits);

				if (kernel.getIRFunction()->nInstructions() > 1) {
					// Execute the kernel!
					kernel(buff);
					nlocalqpucalls++;
				}

				localExpectationValue = buff->getExpectationValueZ();

				// The next iteration will have a different
				// state prep circuit, so toss the current one.
				kernel.getIRFunction()->removeInstruction(0);
			} else {
				localExpectationValue = 1.0;
			}

			// Sum up the expectation values, the Hamiltonian
			// terms coefficient is stored in the first
			// parameter of the Kernels IR Function representation
			expVals.push_back({parameters, std::real(
					boost::get<std::complex<double>>(
							kernel.getIRFunction()->getParameter(0)))
					* localExpectationValue});
		}

	}

	return expVals;
}


}
}
