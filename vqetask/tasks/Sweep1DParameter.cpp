#include "Sweep1DParameter.hpp"
#include "VQEProgram.hpp"


namespace xacc {
namespace vqe {

VQETaskResult Sweep1DParameter::execute(
		Eigen::VectorXd parameters) {
	computeTask = std::make_shared<ComputeEnergyVQETask>(program);
	VQETaskResult result;
	Eigen::VectorXd param(1);
    for (int i = 0; i < parameters.rows(); i++) {
    		param[0] = parameters(i);
    		auto energy = computeTask->execute(param)[0].second;
    		result.push_back({param, energy});
    }
	return result;
}
}
}
