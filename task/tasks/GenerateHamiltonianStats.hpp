#ifndef VQETASKS_GENERATEHAMILTONIANSTATS_HPP_
#define VQETASKS_GENERATEHAMILTONIANSTATS_HPP_

#include "VQETask.hpp"

namespace xacc {
namespace vqe {

class GenerateHamiltonianStats: public VQETask {

public:

	GenerateHamiltonianStats() {}

	GenerateHamiltonianStats(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return the name of this instance.
	 *
	 * @return name The string name
	 */
	virtual const std::string name() const {
		return "vqe-statistics";
	}

	/**
	 * Return the description of this instance
	 * @return description The description of this object.
	 */
	virtual const std::string description() const {
		return "";
	}

};

}
}
#endif
