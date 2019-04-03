#include "ComputeEnergyVQETask.hpp"
#include "IRProvider.hpp"
#include "VQEProgram.hpp"
#include "XACC.hpp"
#include <iomanip>
#include <regex>
#include "xacc_service.hpp"

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

  std::vector<double> vparameters(parameters.data(), parameters.data()+parameters.size());

  // Evaluate our variable parameterized State Prep circuite
  // to produce a state prep circuit with actual rotations
  auto evaluatedStatePrep = statePrep->operator()(vparameters);
  auto optPrep = evaluatedStatePrep->enabledView();

  globalBuffer->addExtraInfo("circuit-depth", optPrep->depth());
  auto qasmStr = optPrep->toString("q");
  qasmStr = std::regex_replace(qasmStr, std::regex("\\n"),"\\\\n");

  globalBuffer->addExtraInfo("ansatz-qasm", qasmStr);

  auto getCoeff = [](Kernel<> &k) -> double {
    return std::real(
        k.getIRFunction()->getParameter(0).as<std::complex<double>>());
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
          std::real(kernel->getParameter(0).as<std::complex<double>>());
      auto expval = b->getExpectationValueZ();
      sum += t * expval;
      b->addExtraInfo("parameters", paramsInfo);
      b->addExtraInfo("kernel", ExtraInfo(kernel->name()));
      b->addExtraInfo("exp-val-z", ExtraInfo(expval));
      b->addExtraInfo("coefficient", ExtraInfo(t));
      globalBuffer->appendChild(kernel->name(), b);
      count++;
    }
    // globalBuffer->appendChild(k)
    xacc::setOption("tnqvm-reset-visitor", "true");
  } else { // NOT TNQVM

    double identityCoeff = 0.0;
    // Create an empty KernelList to be filled
    // with non-trivial kernels
    KernelList<> kernels(qpu);
    // kernels.setBufferPostprocessors(program->getBufferPostprocessors());
    for (auto &k : program->getVQEKernels()) {
      if (k.getIRFunction()->nInstructions() > 0) { // IF NOT IDENTITY TERM
        // If not identity, add the state prep to the circuit
        k.getIRFunction()->insertInstruction(0, optPrep);
        kernels.push_back(k);
      } else { // IF IS IDENTITY TERM
        // if it is identity, add its coeff to the energy sum
        if (rank == 0) {

          identityCoeff = getCoeff(k);

          sum +=identityCoeff;

          auto ibuff = qpu->createBuffer("I", globalBuffer->size());
          ibuff->addExtraInfo("kernel", ExtraInfo("I"));
          ibuff->addExtraInfo("exp-val-z", ExtraInfo(1.0));
          ibuff->addExtraInfo("coefficient", ExtraInfo(identityCoeff));
          ibuff->addExtraInfo("parameters", paramsInfo);
          ibuff->addExtraInfo("ro-fixed-exp-val-z", ExtraInfo(1.0));

          globalBuffer->appendChild("I", ibuff);
        }
      }
    }

    // We can do this in parallel or serially
    if (xacc::optionExists("vqe-use-mpi")) {
      // Allocate some qubits
      auto buf = qpu->createBuffer("tmp", nQubits);
      int myStart = (rank)*kernels.size() / nRanks;
      int myEnd = (rank + 1) * kernels.size() / nRanks;
      for (int i = myStart; i < myEnd; i++) {
        kernels[i](buf);
        totalQpuCalls++;
        sum += getCoeff(kernels[i]) * buf->getExpectationValueZ();
        buf->resetBuffer();
      }

      double result = 0.0;
      int ncalls = 0;
      comm->sumDoubles(sum, result);
      comm->sumInts(totalQpuCalls, ncalls);
      totalQpuCalls = ncalls;
      sum = result;
    } else { // THIS IS NOT TNQVM AND NOT MPI

      // Execute all nontrivial kernels!
      globalBuffer->addExtraInfo("identity-coeff", ExtraInfo(identityCoeff) );
      auto results = kernels.execute(globalBuffer);
      totalQpuCalls += qpu->isRemote() ? 1 : kernels.size();

      // Compute the energy
      for (int i = 0; i < results.size(); ++i) {
        double exp = 0.0;
        if (xacc::optionExists("converge-ro-error") &&
            results[i]->hasExtraInfoKey("ro-fixed-exp-val-z")) {
          exp =
              mpark::get<double>(results[i]->getInformation("ro-fixed-exp-val-z"));
          results[i]->addExtraInfo(
              "exp-val-z", ExtraInfo(results[i]->getExpectationValueZ()));
        } else { // GET THE REGULAR EXP VALUE, NOT RO-FIXED
          exp = results[i]->getExpectationValueZ();
          results[i]->addExtraInfo("exp-val-z", ExtraInfo(exp));
        }

        results[i]->addExtraInfo("parameters", paramsInfo);

        if (results.size() == kernels.size()) {
        auto k = kernels[i];
        auto t = getCoeff(k);
        sum += exp * t;
        results[i]->addExtraInfo("kernel", ExtraInfo(k.getName()));
        results[i]->addExtraInfo("coefficient", ExtraInfo(t));
        globalBuffer->appendChild(k.getName(), results[i]);
        expVals.insert({k.getName(), exp});

        } else {
        auto fname = mpark::get<std::string>(results[i]->getInformation("kernel"));
        globalBuffer->appendChild(fname, results[i]);
        expVals.insert({fname, exp});
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
      [&](ExtraInfo &i) -> bool { return sum < mpark::get<double>(i); });

  if (added) {
    globalBuffer->addExtraInfo("vqe-angles", paramsInfo);
  }

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
