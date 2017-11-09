#ifndef VQETASKS_COMPUTEENERGYVQETASK_HPP_
#define VQETASKS_COMPUTEENERGYVQETASK_HPP_

#include "StatePreparationEvaluator.hpp"
#include "VQETask.hpp"

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

	/**
	 * Return an empty options_description, this is for
	 * subclasses to implement.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		auto desc = std::make_shared<options_description>(
				"Compute Energy VQE Task Options");
		desc->add_options()("vqe-compute-energies-multi-exec",
				"Instead of OpenMP/MPI execution, use XACC multi-execution kernel list.")
				("vqe-compute-persist-buffer-data", value<std::string>(), "Base file name for buffer data.");
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

	int vqeIteration = 0;
	int totalQpuCalls = 0;

};
}
}
#endif
