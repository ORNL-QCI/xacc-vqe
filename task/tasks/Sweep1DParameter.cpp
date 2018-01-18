#include "Sweep1DParameter.hpp"

#include "VQEProgram.hpp"


namespace xacc {
namespace vqe {

VQETaskResult Sweep1DParameter::execute(
		Eigen::VectorXd parameters) {

	computeTask = std::make_shared<ComputeEnergyVQETask>(program);
	VQETaskResult result;
	Eigen::VectorXd param(1);

	if (xacc::optionExists("vqe-restart-index")) {

		int idx = std::stoi(xacc::getOption("vqe-restart-index"));
		for (int i = idx; i < parameters.rows(); i++) {
			param[0] = parameters(i);
			auto energy = computeTask->execute(param).results[0].second;
			result.results.push_back( { param, energy });
		}

	} else {

		for (int i = 0; i < parameters.rows(); i++) {
			param[0] = parameters(i);
			auto energy = computeTask->execute(param).results[0].second;
			result.results.push_back( { param, energy });
		}

	}

	if (xacc::optionExists("vqe-sweep-persist")) {
		std::ofstream out(xacc::getOption("vqe-sweep-persist"));
		for (auto r : result.results) {
			out << r.first(0) << ", " << r.second << "\n";
		}
		out.flush();
		out.close();
	}

	return result;
}
}
}
