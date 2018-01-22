#ifndef IR_VQETASK_HPP_
#define IR_VQETASK_HPP_

#include "OptionsProvider.hpp"
#include <Eigen/Dense>
#include "VQEProgram.hpp"

namespace xacc {
namespace vqe {

class VQETaskResult {
public:
	std::vector<std::pair<Eigen::VectorXd, double>> results;
};

class VQETask : public xacc::Identifiable, public OptionsProvider {

public:

	VQETask() {}

	VQETask(std::shared_ptr<VQEProgram> prog) : program(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters) = 0;

	void setVQEProgram(std::shared_ptr<VQEProgram> p) {
		program = p;
	}

	/**
	 * Return an empty options_description, this is for
	 * subclasses to implement.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		return std::make_shared<options_description>();
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

	virtual ~VQETask() {}

protected:

	std::shared_ptr<VQEProgram> program;
};

}
}

#endif
