#include "ComputeEnergyVQETask.hpp"
#include "XACC.hpp"
#include "VQEProgram.hpp"

namespace xacc {
namespace vqe {

VQETaskResult ComputeEnergyVQETask::execute(
		Eigen::VectorXd parameters) {

	// Local Declarations
	auto comm = program->getCommunicator();
	VQETaskResult taskResult;
	double sum = 0.0, localExpectationValue = 0.0;
	int rank = comm->rank(), nlocalqpucalls = 0;
	int nRanks = comm->size();
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

	// Separate our kernels into Identity kernels and 
	// non trivial kernels to be executed
	KernelList<> identityKernels, kernels(qpu);
	std::map<std::string, double> nameToCoeff;
	std::map<std::string, std::string> nameToTag;
	for (auto& k : program->getVQEKernels()) {
		if (k.getIRFunction()->nInstructions() > 0) {
			k.getIRFunction()->insertInstruction(0,evaluatedStatePrep);
			kernels.push_back(k);
			nameToTag.insert({k.getName(), k.getIRFunction()->getTag()});
			double coeff = std::real(boost::get<std::complex<double>>(k.getIRFunction()->getParameter(0)));
			nameToCoeff.insert({k.getName(), coeff});
		} else {
			identityKernels.push_back(k);
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
			if(!isReadoutErrorKernel(kernels[i].getIRFunction()->getTag()))
				sum += nameToCoeff[kernels[i].getName()] * buf->getExpectationValueZ();
			buf->resetBuffer();
		}

		if (rank == 0) 
			for (auto& idK : identityKernels) sum += std::real(boost::get<std::complex<double>>(idK.getIRFunction()->getParameter(0)));

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
		int i = 0;
		for(auto& kv : nameToCoeff) {
			std::cout << "Computing " << kv.first << ", " << kv.second << ", " << results[i]->getExpectationValueZ() << "\n";
			if(!isReadoutErrorKernel(nameToTag[kv.first])) {
				sum += results[i]->getExpectationValueZ() * kv.second;
				++i;
			}
		}
		for(auto& idK : identityKernels) sum += std::real(boost::get<std::complex<double>>(idK.getIRFunction()->getParameter(0)));
		for(auto& k : kernels) k.getIRFunction()->removeInstruction(0);	
	}

	std::stringstream ss;
	ss << std::setprecision(10) << sum << " at (" << parameters.transpose() << ")";
	if (rank == 0) {
		xacc::info("Iteration " + std::to_string(vqeIteration) + ", Computed VQE Energy = " + ss.str());
	}

	vqeIteration++;
	taskResult.results.push_back({parameters, sum});
	taskResult.energy = sum;
	taskResult.angles = parameters;
	taskResult.nQpuCalls = totalQpuCalls;

//	std::stringstream outputString;
//	if (persist) {
//		if (!boost::filesystem::exists(xacc::getOption("vqe-persist-data"))) {
//			for (int i = 0; i < program->getNParameters(); i++) {
//				outputString << "theta" << i << ", ";
//			}
//
//			for (int i = 0; i < kernelNames.size(); i++) {
//				outputString << (i == 0 ? "" : ", ") << kernelNames[i];
//			}
//
//			outputString << ", E\n";
//		}
//		for (int i = 0; i < parameters.rows(); i++) {
//			outputString << parameters(i) << ", ";
//		}
//	}

//	double maxExecTime = 0.0;

	return taskResult;
}


}
}
