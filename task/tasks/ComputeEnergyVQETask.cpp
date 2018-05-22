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
	std::map<std::string, double> expVals;
	bool multiExec = false;
	bool persist = xacc::optionExists("vqe-persist-data");

	auto statePrep = program->getStatePreparationCircuit();
	auto nQubits = program->getNQubits();
	auto qpu = program->getAccelerator();

	// Evaluate our variable parameterized State Prep circuite
	// to produce a state prep circuit with actual rotations
	auto evaluatedStatePrep = StatePreparationEvaluator::evaluateCircuit(
			statePrep, program->getNParameters(), parameters);

	auto isReadoutErrorKernel = [](const std::string& tag) -> bool 
		{ return boost::contains(tag,"readout-error");};

	auto getCoeff = [](Kernel<>& k) ->double {
		return std::real(boost::get<std::complex<double>>(k.getIRFunction()->getParameter(0)));
	};

	// Separate our kernels into Identity kernels and 
	// non trivial kernels to be executed
	KernelList<> kernels(qpu);

	double idSum = 0.0;
	for (auto& k : program->getVQEKernels()) {
		if (k.getIRFunction()->nInstructions() > 0) {
			k.getIRFunction()->insertInstruction(0,evaluatedStatePrep);
			kernels.push_back(k);
		} else {
			if (rank == 0) sum += getCoeff(k);
		}
	}

	// Allocate some qubits
	auto buf = qpu->createBuffer("tmp",nQubits);

	if (xacc::optionExists("vqe-use-mpi")) {
		int myStart = (rank) * kernels.size() / nRanks;
		int myEnd = (rank + 1) * kernels.size() / nRanks;
		for (int i = myStart; i < myEnd; i++) {
			kernels[i](buf);
			totalQpuCalls++;
			if(!isReadoutErrorKernel(kernels[i].getIRFunction()->getTag())) {
				sum += getCoeff(kernels[i]) 
					* buf->getExpectationValueZ();
			}
			buf->resetBuffer();
		}

		double result = 0.0;
		int ncalls = 0;
		comm->sumDoubles(sum, result);
		comm->sumInts(totalQpuCalls, ncalls);
		totalQpuCalls = ncalls;
		sum = result;
	} else {
		// Execute!
		auto results = kernels.execute(buf);
		totalQpuCalls += qpu->isRemote() ? 1 : kernels.size();

		// Compute the energy
		for(int i = 0; i < results.size(); ++i) {
			auto k = kernels[i];
			auto exp = results[i]->getExpectationValueZ();
			if(!isReadoutErrorKernel(k.getIRFunction()->getTag())) {
				sum += exp * 
					getCoeff(k);
			}
			expVals.insert({k.getName(), exp});
		}
		for(auto& k : kernels) k.getIRFunction()->removeInstruction(0);	
	}

	std::stringstream ss;
	ss << std::setprecision(10) << sum << " at (" << parameters.transpose() << ")";
	if (rank == 0) {
		xacc::info("Iteration " + std::to_string(vqeIteration) + ", Computed VQE Energy = " + ss.str());
	}

	vqeIteration++;

	if (persist) {
		VQETaskResult taskResult(xacc::getOption("vqe-persist-data"));
		taskResult.energy = sum;
		taskResult.angles = parameters;
		taskResult.nQpuCalls = totalQpuCalls;
		taskResult.expVals = expVals;
		taskResult.persist();
		return taskResult;
	} else {
		VQETaskResult taskResult;
		taskResult.energy = sum;
		taskResult.angles = parameters;
		taskResult.nQpuCalls = totalQpuCalls;
		taskResult.expVals = expVals;
		return taskResult;
	}
}


}
}
