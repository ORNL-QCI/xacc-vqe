#ifndef TASK_TASKS_ARPACK_SLEPCGROUNDSTATEENERGYCALCULATOR_HPP_
#define TASK_TASKS_ARPACK_SLEPCGROUNDSTATEENERGYCALCULATOR_HPP_

#include "BruteForceComputeGroundStateEnergy.hpp"

namespace xacc {
namespace vqe {

class SlepcGroundStateEnergyCalculator: public GroundStateEnergyCalculator {
	virtual double computeGroundStateEnergy(CompositeSpinInstruction& inst,
			const int nQubits);

	virtual const std::string name() const {
		return "vqe-slepc";
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
