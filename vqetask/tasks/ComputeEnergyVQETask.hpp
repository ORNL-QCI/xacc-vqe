#ifndef VQETASKS_COMPUTEENERGYVQETASK_HPP_
#define VQETASKS_COMPUTEENERGYVQETASK_HPP_

#include "VQETask.hpp"
#include "StatePreparationEvaluator.hpp"

namespace xacc {
namespace vqe {

class ComputeEnergyVQETask: public VQETask {

public:

	ComputeEnergyVQETask() {}

	ComputeEnergyVQETask(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return the name of this instance.
	 *
	 * @return name The string name
	 */
	virtual const std::string name() const {
		return "compute-energy";
	}

	/**
	 * Return the description of this instance
	 * @return description The description of this object.
	 */
	virtual const std::string description() const {
		return "";
	}

	int vqeIteration = 0;
	int totalQpuCalls = 0;

};
}
}
#endif
