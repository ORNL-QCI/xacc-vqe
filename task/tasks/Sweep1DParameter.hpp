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

	/**
	 * Return an empty options_description, this is for
	 * subclasses to implement.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		auto desc = std::make_shared<options_description>(
				"Sweep 1D Energies");
		desc->add_options()("vqe-sweep-persist",value<std::string>(),
				"Persist the results to the provided CSV file.")
				("vqe-restart-index", value<std::string>(), "");
		return desc;
	}

	/**
	 * Given user-input command line options, perform
	 * some operation. Returns true if runtime should exit,
	 * false otherwise.
	 *
	 * @param map The mapping of options to values
	 * @return exit True if exit, false otherwise
	 */
	virtual bool handleOptions(variables_map& map) {
		return false;
	}

	std::shared_ptr<ComputeEnergyVQETask> computeTask;

};
}
}
#endif
