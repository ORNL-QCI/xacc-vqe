#ifndef VQETASKS_VQEMINIMIZETASK_HPP_
#define VQETASKS_VQEMINIMIZETASK_HPP_

#include "VQETask.hpp"
#include "ComputeEnergyVQETask.hpp"

#include "solver/neldermeadsolver.h"
#include "solver/conjugatedgradientdescentsolver.h"
#include "solver/gradientdescentsolver.h"

namespace xacc {
namespace vqe {

class VQEMinimizeTask: public VQETask, public cppoptlib::Problem<double> {

public:

	VQEMinimizeTask() {}

	VQEMinimizeTask(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return the name of this instance.
	 *
	 * @return name The string name
	 */
	virtual const std::string name() const {
		return "vqe";
	}

	/**
	 * Return the description of this instance
	 * @return description The description of this object.
	 */
	virtual const std::string description() const {
		return "";
	}

	double value(const Eigen::VectorXd& x);

	/**
	 * This internal structure describes the criteria for
	 * convergence in this VQE execution.
	 */
	class VQECriteria: public cppoptlib::Criteria<double> {
	public:
		static VQECriteria defaults() {
			VQECriteria d;
			d.iterations = 1000;
			d.xDelta = 0;
			d.fDelta = 1e-6;
			d.gradNorm = 1e-4;
			d.condition = 0;
			return d;
		}
	};

	/**
	 * This static method returns the user-specified
	 * convergence criteria.
	 *
	 * @return criteria The convergence criteria.
	 */
	static VQECriteria getConvergenceCriteria() {
		auto criteria = VQECriteria::defaults();

		if (xacc::optionExists("vqe-energy-delta")) {
			criteria.fDelta = std::stod(xacc::getOption("vqe-energy-delta"));
		}

		if (xacc::optionExists("vqe-iterations")) {
			criteria.iterations = std::stoi(xacc::getOption("vqe-iterations"));
		}

		return criteria;

	}

protected:

	double currentEnergy = 0.0;

	std::shared_ptr<ComputeEnergyVQETask> computeTask;

};
}
}
#endif
