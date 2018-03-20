#ifndef TASK_TASKS_PSOVQEBACKEND_HPP_
#define TASK_TASKS_PSOVQEBACKEND_HPP_

#include "VQEMinimizeTask.hpp"

namespace xacc {
namespace vqe {

class PsoVQEBackend: public VQEBackend, public OptionsProvider {

protected:

	std::shared_ptr<ComputeEnergyVQETask> computeTask;

public:

	virtual const VQETaskResult minimize(Eigen::VectorXd parameters);

	virtual const std::string name() const {
		return "pso";
	}

	/**
	 * Return the description of this instance
	 * @return description The description of this object.
	 */
	virtual const std::string description() const {
		return "";
	}

	/**
	 * Return an empty options_description, this is for
	 * subclasses to implement.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		auto desc = std::make_shared<options_description>(
				"Pso Options");
		return desc;
	}

	virtual bool handleOptions(variables_map& map) {
		return false;
	}

};

}
}
#endif
