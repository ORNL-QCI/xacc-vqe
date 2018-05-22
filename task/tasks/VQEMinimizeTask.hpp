#ifndef VQETASKS_VQEMINIMIZETASK_HPP_
#define VQETASKS_VQEMINIMIZETASK_HPP_

#include "ComputeEnergyVQETask.hpp"
#include "VQETask.hpp"
#include "solver/neldermeadsolver.h"
#include "solver/conjugatedgradientdescentsolver.h"
#include "solver/gradientdescentsolver.h"
#include "OptionsProvider.hpp"

namespace xacc {
namespace vqe {

class VQEBackend : public Identifiable {
protected:

	std::shared_ptr<VQEProgram> program;
public:
	virtual const VQETaskResult minimize(Eigen::VectorXd parameters) = 0;
	virtual void setProgram(std::shared_ptr<VQEProgram> p) { program = p;}
	virtual ~VQEBackend(){}
};

class CppOptVQEBackend : public VQEBackend, public cppoptlib::Problem<double> {

protected:

	double currentEnergy = 0.0;

	std::shared_ptr<ComputeEnergyVQETask> computeTask;

public:

	CppOptVQEBackend() {
	}

	virtual const VQETaskResult minimize(Eigen::VectorXd parameters) {
		computeTask = std::make_shared<ComputeEnergyVQETask>(program);
		cppoptlib::NelderMeadSolver<CppOptVQEBackend> solver;
		solver.setStopCriteria(CppOptVQEBackend::getConvergenceCriteria());
		solver.minimize(*this, parameters);
		VQETaskResult result;
		result.angles = parameters;
		result.energy = currentEnergy;
		result.nQpuCalls = computeTask->totalQpuCalls;
		result.vqeIterations = computeTask->vqeIteration;
		return result;
	}

	double value(const Eigen::VectorXd& x) {
		currentEnergy = computeTask->execute(x).energy;
		return currentEnergy;
	}

	virtual const std::string name() const {
		return "cppopt";
	}

	/**
	 * Return the description of this instance
	 * @return description The description of this object.
	 */
	virtual const std::string description() const {
		return "";
	}


	/**
	 * This internal structure describes the criteria for
	 * convergence in this VQE execution.
	 */
	class VQECriteria: public cppoptlib::Criteria<double> {
	public:
		static VQECriteria defaults() {
			VQECriteria d;
			d.iterations = 10000;
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
};

class VQEMinimizeTask: public VQETask {

public:

	VQEMinimizeTask() {}

	VQEMinimizeTask(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return an empty options_description, this is for
	 * subclasses to implement.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		auto desc = std::make_shared<options_description>(
				"VQE Task Options");
		desc->add_options()("vqe-backend", value<std::string>(),
							"The backend to use to compute the min energy via VQE");
		return desc;
	}

	virtual bool handleOptions(variables_map& map) {
		return false;
	}
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

};
}
}
#endif
