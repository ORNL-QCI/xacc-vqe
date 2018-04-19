#include "ComputeEnergyVQETask.hpp"
#include "XACC.hpp"
#include "VQEProgram.hpp"

namespace xacc {
namespace vqe {

VQETaskResult ComputeEnergyVQETask::execute(
		Eigen::VectorXd parameters) {

	// Local Declarations
	auto comm = program->getCommunicator();
	double sum = 0.0, localExpectationValue = 0.0;
	int rank = comm->rank(), nlocalqpucalls = 0;
	int nRanks = comm->size();
	bool multiExec = false;
	bool persist = xacc::optionExists("vqe-persist-data");

	auto statePrep = program->getStatePreparationCircuit();
	auto nQubits = program->getNQubits();
	auto qpu = program->getAccelerator();
	int nReadoutKernels = 0;

	auto pauli = program->getPauliOperator();

	// Evaluate our variable parameterized State Prep circuite
	// to produce a state prep circuit with actual rotations
	auto evaluatedStatePrep = StatePreparationEvaluator::evaluateCircuit(
			statePrep, program->getNParameters(), parameters);

	auto kernels = program->getVQEKernels();

	VQETaskResult taskResult;

	std::vector<std::string> kernelNames;
	for (auto& k : kernels) {
		if (k.getIRFunction()->getTag() == "readout-error") {
			nReadoutKernels++;
		} else if (k.getIRFunction()->nInstructions() > 0){
			kernelNames.push_back(k.getIRFunction()->name());
		}
	}

	if (xacc::optionExists("vqe-compute-energies-multi-exec")
			|| qpu->name() == "ibm") {
		multiExec = true;
	}

	std::stringstream outputString;

	if (persist) {

		if (!boost::filesystem::exists(xacc::getOption("vqe-persist-data"))) {
			for (int i = 0; i < program->getNParameters(); i++) {
				outputString << "theta" << i << ", ";
			}

			for (int i = 0; i < kernelNames.size(); i++) {
				outputString << (i == 0 ? "" : ", ") << kernelNames[i];
			}

			outputString << ", E\n";
		}
		for (int i = 0; i < parameters.rows(); i++) {
			outputString << parameters(i) << ", ";
		}
	}

	double maxExecTime = 0.0;

	if (multiExec) {

		xacc::info("Computing Energy with XACC Kernel Multi-Exec.");
		std::vector<double> coeffs;
		std::vector<Kernel<>> identityKernels;

		for (int i = 0; i < kernels.size(); i++) {
			auto k = kernels[i];
			// If not identity
			if (k.getIRFunction()->nInstructions() > 0
					&& !boost::contains(k.getIRFunction()->getTag(),
							"readout-error")) {
				k.getIRFunction()->insertInstruction(0, evaluatedStatePrep);
				coeffs.push_back(
						std::real(
								boost::get<std::complex<double>>(
										k.getIRFunction()->getParameter(0))));
			} else if(!boost::contains(k.getIRFunction()->getTag(),
							"readout-error")) {
				identityKernels.push_back(k);
			}
		}

		auto buff = qpu->createBuffer("qreg", nQubits);


		// FIXME, Goal here is to run AcceleratorBufferPostprocessors in KernelList.execute...
		auto tmpBuffers = kernels.execute(buff);
		nlocalqpucalls++;

		maxExecTime += qpu->getExecutionTime();

		execTime += maxExecTime;

		int counter = 0;
		std::vector<double> expVals;
		for (int i = 0; i < coeffs.size(); i++) {
			localExpectationValue = tmpBuffers[i]->getExpectationValueZ();
			xacc::info(
					"Fixed Expectation for Kernel " + std::to_string(counter) + " = "
							+ std::to_string(localExpectationValue));
			expVals.push_back(localExpectationValue);
			sum += coeffs[counter] * localExpectationValue;
			counter++;
		}

		for (auto k : identityKernels) {
			sum += std::real(
					boost::get<std::complex<double>>(
							k.getIRFunction()->getParameter(0)));
		}

		if (persist) {

			for (int i = 0; i < expVals.size(); i++) {
				outputString << (i==0 ? "" : ", ") << expVals[i];
			}

			outputString << ", " << sum << "\n";

		}
		for (int i = 0; i < kernels.size(); i++) {
			if (kernels[i].getIRFunction()->nInstructions() > 0) {
				kernels[i].getIRFunction()->removeInstruction(0);
			}
		}
	} else {

		std::map<std::string, double> expVals;

		for (auto kn : kernelNames) {
			expVals.insert({kn, 0.0});
		}

		double localExecTime = 0.0;

		// Execute the kernels on the appropriate QPU
		// in parallel MPI, get local rank's Kernels to execute
		int myStart = (rank) * kernels.size() / nRanks;
		int myEnd = (rank + 1) * kernels.size() / nRanks;

		for (int i = myStart; i < myEnd; i++) {
			
			double lexpval = 1.0;
			
			// Get the ith Kernel
			auto kernel = kernels[i];

			if (!boost::contains(kernel.getIRFunction()->getTag(),
					"readout-error")) {
				auto t = std::real(
						boost::get<std::complex<double>>(
								kernel.getIRFunction()->getParameter(0)));

				// If not an identity kernel...
				if (kernel.getIRFunction()->nInstructions() > 0) {

					// Insert the state preparation circuit IR
					// at location 0 in this Kernels IR instructions.
					kernel.getIRFunction()->insertInstruction(0,
							evaluatedStatePrep);

					// Create a temporary buffer of qubits
					auto buff = qpu->createBuffer("qreg", nQubits);

					// Execute the kernel!
					kernel(buff);
					nlocalqpucalls++;

					localExecTime += qpu->getExecutionTime();

					lexpval = buff->getExpectationValueZ();

					expVals[kernel.getIRFunction()->name()] = lexpval;

					// The next iteration will have a different
					// state prep circuit, so toss the current one.
					kernel.getIRFunction()->removeInstruction(0);
				}

				// Sum up the expectation values, the Hamiltonian
				// terms coefficient is stored in the first
				// parameter of the Kernels IR Function representation
				sum += t * lexpval;
			}
		}

		// Find the max localExecTime across all ranks
		comm->maxDouble(localExecTime, maxExecTime);
		execTime += maxExecTime;

		if (persist) {

			int counter = 0;
			for (auto& kn : kernelNames) {
				outputString << (counter == 0 ? "" : ", ") << expVals[kn];
				counter++;
			}

			outputString << ", " << sum << "\n";
		}
	}

	double result = 0.0;
	int totalqpucalls = 0;
	comm->sumDoubles(sum, result);
	comm->sumInts(nlocalqpucalls, totalqpucalls);

	double currentEnergy = result;
	totalQpuCalls += totalqpucalls;

	std::stringstream ss;
	ss << std::setprecision(10) << currentEnergy << " at (" << parameters.transpose() << "), exec_time = " << maxExecTime << " (s)";

	if (rank == 0) {
		xacc::info("Iteration " + std::to_string(vqeIteration) + ", Computed VQE Energy = " + ss.str());
	}

	vqeIteration++;

	taskResult.results.push_back({parameters, currentEnergy});

	taskResult.energy = currentEnergy;
	taskResult.angles = parameters;
	taskResult.nQpuCalls = totalQpuCalls;
	taskResult.execTime = execTime;

	if (persist) {
		std::ofstream out(xacc::getOption("vqe-persist-data"), std::ios_base::app);
		out << outputString.str();
	}

	return taskResult;
}


}
}
