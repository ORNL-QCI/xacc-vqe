#include "VQEMinimizeTask.hpp"

#include "VQEProgram.hpp"


namespace xacc {
namespace vqe {


VQETaskResult VQEMinimizeTask::execute(
		Eigen::VectorXd parameters) {

	auto serviceReg = ServiceRegistry::instance();

	std::shared_ptr<VQEBackend> backend;
	if (xacc::optionExists("vqe-backend")) {
		backend = serviceReg->getService<VQEBackend>(xacc::getOption("vqe-backend"));
	} else {
		backend = std::make_shared<CppOptVQEBackend>();
	}

	backend->setProgram(program);
	auto result = backend->minimize(parameters);

	std::stringstream ss;
	ss << result.nQpuCalls << " total QPU calls over " << result.vqeIterations << " VQE iterations.";
	xacc::info("");
	xacc::info(ss.str());
	return result;
}

}
}
