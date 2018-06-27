#include "VQEMinimizeTask.hpp"

#include "VQEProgram.hpp"


namespace xacc {
namespace vqe {


VQETaskResult VQEMinimizeTask::execute(
		Eigen::VectorXd parameters) {

	std::shared_ptr<VQEBackend> backend;
	if (xacc::optionExists("vqe-backend")) {
		backend = xacc::getService<VQEBackend>(xacc::getOption("vqe-backend"));
	} else {
		backend = std::make_shared<CppOptVQEBackend>();
	}

	backend->setProgram(program);
	auto f = program->getStatePreparationCircuit();
	auto result = backend->minimize(parameters);
	result.ansatzQASM = f->toString("q");
	std::stringstream ss;
	ss << result.nQpuCalls << " total QPU calls over " << result.vqeIterations << " VQE iterations.";
	xacc::info("");
	xacc::info(ss.str());
	return result;
}

}
}
