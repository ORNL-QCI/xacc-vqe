#ifndef VQETASKS_SWEEP1DPARAMETER_HPP_
#define VQETASKS_SWEEP1DPARAMETER_HPP_

#include "../../task/tasks/ComputeEnergyVQETask.hpp"
#include "../../task/VQETask.hpp"

namespace xacc {
namespace vqe {

class Sweep1DParameter: public VQETask {

public:

	Sweep1DParameter() {}

	Sweep1DParameter(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return the name of this instance.
	 *
	 * @return name The string name
	 */
	virtual const std::string name() const {
		return "sweep-1d";
	}

	/**
	 * Return the description of this instance
	 * @return description The description of this object.
	 */
	virtual const std::string description() const {
		return "";
	}

	std::shared_ptr<ComputeEnergyVQETask> computeTask;

};
}
}
#endif
