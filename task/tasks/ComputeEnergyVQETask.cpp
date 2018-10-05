#include "ComputeEnergyVQETask.hpp"
#include "IRProvider.hpp"
#include "VQEProgram.hpp"
#include "XACC.hpp"
#include <boost/algorithm/string.hpp>

namespace xacc {
namespace vqe {

VQETaskResult ComputeEnergyVQETask::execute(Eigen::VectorXd parameters) {

  // Local Declarations
  auto comm = program->getCommunicator();
  double sum = 0.0;
  int rank = comm->rank(), nlocalqpucalls = 0;
  int nRanks = comm->size();
  std::map<std::string, double> expVals, readoutProbs;
  bool persist = xacc::optionExists("vqe-persist-data");

  auto globalBuffer = program->getGlobalBuffer();
  std::vector<double> paramsVec(parameters.size());
  for (int i = 0; i < parameters.size(); i++) {
    paramsVec[i] = parameters(i);
  }
  ExtraInfo paramsInfo(paramsVec);

  // Get info about the problem
  auto statePrep = program->getStatePreparationCircuit();
  auto nQubits = program->getNQubits();
  auto qpu = program->getAccelerator();

  // Evaluate our variable parameterized State Prep circuite
  // to produce a state prep circuit with actual rotations
  auto evaluatedStatePrep = statePrep->operator()(parameters);
  auto optPrep = evaluatedStatePrep->enabledView();

  globalBuffer->addExtraInfo("circuit-depth",optPrep->depth());
  auto qasmStr = optPrep->toString("q");
  boost::replace_all(qasmStr, "\\n","\\\\n");
  globalBuffer->addExtraInfo("ansatz-qasm", qasmStr);
  
//   xacc::info("StatePrep:\n" + optPrep->toString("q"));
  // Utility functions for readability
  auto isReadoutErrorKernel = [](const std::string &tag) -> bool {
    return boost::contains(tag, "readout-error");
  };
  auto getCoeff = [](Kernel<> &k) -> double {
    return std::real(
        boost::get<std::complex<double>>(k.getIRFunction()->getParameter(0)));
  };

  if (qpu->name() == "tnqvm" && !xacc::optionExists("vqe-use-mpi")) {
    xacc::setOption("run-and-measure", "");
    std::vector<std::shared_ptr<Function>> ks;
    ks.push_back(optPrep);
    for (auto &k : program->getVQEKernels()) {
      ks.push_back(k.getIRFunction());
    }
    auto buffer = qpu->createBuffer("q", nQubits);
    auto tmpBuffers = qpu->execute(buffer, ks);
    ks.erase(ks.begin());
    int count = 0;
    for (auto &b : tmpBuffers) {
      auto kernel = ks[count];
      auto t =
          std::real(boost::get<std::complex<double>>(kernel->getParameter(0)));
      auto expval = b->getExpectationValueZ();
      sum += t * expval;
      b->addExtraInfo("parameters", paramsInfo);
      b->addExtraInfo("kernel", ExtraInfo(kernel->name()));
      b->addExtraInfo("exp-val-z", ExtraInfo(expval));
      
      globalBuffer->appendChild(kernel->name(), b);
      count++;
    }
    // globalBuffer->appendChild(k)
    xacc::setOption("tnqvm-reset-visitor", "true");
  } else {
    // Create an empty KernelList to be filled
    // with non-trivial kernels
    KernelList<> kernels(qpu);
    // kernels.setBufferPostprocessors(program->getBufferPostprocessors());
    for (auto &k : program->getVQEKernels()) {
      if (k.getIRFunction()->nInstructions() > 0) {
        // If not identity, add the state prep to the circuit
        if (k.getIRFunction()->getTag() != "readout-error") {
          k.getIRFunction()->insertInstruction(0, optPrep);
        }
        kernels.push_back(k);
      } else {
        // if it is identity, add its coeff to the energy sum
        if (rank == 0)
          sum += getCoeff(k);
      }
    }

    // Allocate some qubits
    auto buf = qpu->createBuffer("tmp", nQubits);

    // We can do this in parallel or serially
    if (xacc::optionExists("vqe-use-mpi")) {
      int myStart = (rank)*kernels.size() / nRanks;
      int myEnd = (rank + 1) * kernels.size() / nRanks;
      for (int i = myStart; i < myEnd; i++) {
        kernels[i](buf);
        totalQpuCalls++;
        if (!isReadoutErrorKernel(kernels[i].getIRFunction()->getTag())) {
          sum += getCoeff(kernels[i]) * buf->getExpectationValueZ();
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
      // Execute all nontrivial kernels!
      auto results = kernels.execute(buf);
      totalQpuCalls += qpu->isRemote() ? 1 : kernels.size();

      // Compute the energy
      for (int i = 0; i < results.size(); ++i) {
        auto k = kernels[i];
        auto exp = results[i]->getExpectationValueZ();
        if (!isReadoutErrorKernel(k.getIRFunction()->getTag())) {
          sum += exp * getCoeff(k);
          results[i]->addExtraInfo("parameters", paramsInfo);
          results[i]->addExtraInfo("kernel", ExtraInfo(k.getName()));
          results[i]->addExtraInfo("exp-val-z", ExtraInfo(exp));
          globalBuffer->appendChild(k.getName(), results[i]);
        }
        if (k.getIRFunction()->getTag() != "readout-error") {
          expVals.insert({k.getName(), exp});
        } else if (xacc::optionExists("correct-readout-errors")) {
          auto name = k.getName(); // Qubit_{0,1}
          std::vector<std::string> split;
          boost::split(split, name, boost::is_any_of("_"));
          auto qbit = std::stoi(split[0]);
          auto zero_one = std::stoi(split[1]);
          std::string bitstr = "";
          for (int i = 0; i < nQubits; i++)
            bitstr += "0";
          if (zero_one == 0)
            bitstr[nQubits - qbit - 1] = '1';
          auto prob = results[i]->computeMeasurementProbability(bitstr);
          readoutProbs.insert({"p_" + name, std::isnan(prob) ? 0.0 : prob});
        }
      }

      // Clean up by removing the state prep
      // from the measurement kernels
      for (auto &k : kernels)
        k.getIRFunction()->removeInstruction(0);
    }
  }

  std::stringstream ss;
  ss << std::setprecision(10) << sum << " at (" << parameters.transpose()
     << ")";
  if (rank == 0) {
    xacc::info("Iteration " + std::to_string(vqeIteration) +
               ", Computed VQE Energy = " + ss.str());
  }

  vqeIteration++;

  auto added = globalBuffer->addExtraInfo(
      "vqe-energy", ExtraInfo(sum),
      [&](ExtraInfo &i) -> bool { return sum < boost::get<double>(i); });
  if (added)
    globalBuffer->addExtraInfo("vqe-angles", paramsInfo);
  globalBuffer->addExtraInfo("vqe-nQPU-calls", ExtraInfo(totalQpuCalls));

  // See if the user requested data persisitence
  if (persist) {
    VQETaskResult taskResult(xacc::getOption("vqe-persist-data"));
    taskResult.energy = sum;
    taskResult.angles = parameters;
    taskResult.nQpuCalls = totalQpuCalls;
    taskResult.expVals = expVals;
    taskResult.readoutErrorProbabilities = readoutProbs;
    taskResult.persist();
    return taskResult;
  } else {
    VQETaskResult taskResult;
    taskResult.energy = sum;
    taskResult.angles = parameters;
    taskResult.nQpuCalls = totalQpuCalls;
    taskResult.expVals = expVals;
    taskResult.readoutErrorProbabilities = readoutProbs;
    return taskResult;
  }
}

} // namespace vqe
} // namespace xacc
