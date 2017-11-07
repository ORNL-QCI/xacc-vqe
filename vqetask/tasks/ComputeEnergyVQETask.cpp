#include "ComputeEnergyVQETask.hpp"
#include "VQEProgram.hpp"

namespace xacc {
namespace vqe {

VQETaskResult ComputeEnergyVQETask::execute(
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
//	XACCInfo("HELLO WORLD:\n" + evaluatedStatePrep->toString("qreg"));

	auto kernels = program->getVQEKernels();

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

		int counter = 0;
		for (auto b : tmpBuffers) {
			localExpectationValue = b->getExpectationValueZ();
			sum += coeffs[counter] * localExpectationValue;
			counter++;
		}

		for (auto k : identityKernels) {
			sum += std::real(
					boost::get<std::complex<double>>(
							k.getIRFunction()->getParameter(0)));
		}

		for (int i = 0; i < kernels.size(); i++) {
			if (kernels[i].getIRFunction()->nInstructions() > 0) {
				kernels[i].getIRFunction()->removeInstruction(0);
			}
		}
	} else {

		// Execute the kernels on the appropriate QPU
		// in parallel using OpenMP threads per
		// every MPI rank.
		int myStart = (rank) * kernels.size() / nRanks;
		int myEnd = (rank + 1) * kernels.size() / nRanks;
#pragma omp parallel for reduction (+:sum, nlocalqpucalls)
		for (int i = myStart; i < myEnd; i++) {

			// Get the ith Kernel
			auto kernel = kernels[i];

			// Insert the state preparation circuit IR
			// at location 0 in this Kernels IR instructions.
			kernel.getIRFunction()->insertInstruction(0, evaluatedStatePrep);

			// Create a temporary buffer of qubits
			auto buff = qpu->createBuffer("qreg", nQubits);

			if (kernel.getIRFunction()->nInstructions() > 1) {
				// Execute the kernel!
				kernel(buff);
				nlocalqpucalls++;
			}

			// Get Expectation value. The second parameter of
			// the Kernel's IR function stores whether or not
			// this kernel is the Identity.
			if (kernel.getIRFunction()->nInstructions() == 0) {
				localExpectationValue = 1.0;
			} else {
				localExpectationValue = buff->getExpectationValueZ();
			}

			// Sum up the expectation values, the Hamiltonian
			// terms coefficient is stored in the first
			// parameter of the Kernels IR Function representation
			sum += std::real(
					boost::get<std::complex<double>>(
							kernel.getIRFunction()->getParameter(0)))
					* localExpectationValue;

			// The next iteration will have a different
			// state prep circuit, so toss the current one.
			kernel.getIRFunction()->removeInstruction(0);
		}

	}

	// Set the energy.
	double currentEnergy = sum;

	double result = 0.0;
	int totalqpucalls = 0;
	boost::mpi::all_reduce(comm, currentEnergy, result, std::plus<double>());
	boost::mpi::all_reduce(comm, nlocalqpucalls, totalqpucalls, std::plus<double>());
	currentEnergy = result;

	totalQpuCalls += nlocalqpucalls;

	std::stringstream ss;
	ss << currentEnergy << " at (" << parameters.transpose() << ")";

	if (rank == 0) {
		XACCInfo("Iteration " + std::to_string(vqeIteration) + ", Computed VQE Energy = " + ss.str());
//		XACCInfo("\tTotal QPU calls = " + std::to_string(totalQpuCalls));
	}

	vqeIteration++;
	return std::vector<std::pair<Eigen::VectorXd, double>> { { parameters,
			currentEnergy } };
}


}
}
