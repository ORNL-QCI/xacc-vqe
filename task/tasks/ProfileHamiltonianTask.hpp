#ifndef VQETASKS_PROFILEHAMILTONIANTASK_HPP_
#define VQETASKS_PROFILEHAMILTONIANTASK_HPP_

#include "VQETask.hpp"

namespace xacc {
namespace vqe {

class ProfileHamiltonianTask : public VQETask {

public:

	ProfileHamiltonianTask() {}

	ProfileHamiltonianTask(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return the name of this instance.
	 *
	 * @return name The string name
	 */
	virtual const std::string name() const {
		return "vqe-profile";
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
				"Profile Hamiltonian");
		desc->add_options()("vqe-profile-name", value<std::string>(), "The name of the file to save.");
		return desc;
	}

	virtual bool handleOptions(variables_map& map) {
		return false;
	}

};

}
}
#endif
