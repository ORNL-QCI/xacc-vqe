#ifndef TASK_TASKS_SLEPCDIAGONALIZEBACKEND_HPP_
#define TASK_TASKS_SLEPCDIAGONALIZEBACKEND_HPP_

#include "DiagonalizeTask.hpp"

namespace xacc {
namespace vqe {

class SlepcDiagonalizeBackend: public DiagonalizeBackend {
	virtual double diagonalize(std::shared_ptr<VQEProgram> prog);

	virtual const std::string name() const {
		return "slepc";
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
