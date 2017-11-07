#ifndef IR_VQETASK_HPP_
#define IR_VQETASK_HPP_

#include "VQEProgram.hpp"
#include <Eigen/Dense>

namespace xacc {
namespace vqe {

using VQETaskResult = std::vector<std::pair<Eigen::VectorXd, double>>;

class VQETask : public xacc::Identifiable {

public:

	VQETask() {}

	VQETask(std::shared_ptr<VQEProgram> prog) : program(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters) = 0;

	void setVQEProgram(std::shared_ptr<VQEProgram> p) {
		program = p;
	}

	virtual ~VQETask() {}

protected:

	std::shared_ptr<VQEProgram> program;
};

}
}

#endif
