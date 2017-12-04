#ifndef TASK_TASKS_ARPACK_ARPACKPPGROUNDSTATEENERGYCALCULATOR_HPP_
#define TASK_TASKS_ARPACK_ARPACKPPGROUNDSTATEENERGYCALCULATOR_HPP_

#include "BruteForceComputeGroundStateEnergy.hpp"

namespace xacc {
namespace vqe {

class ArpackppGroundStateEnergyCalculator: public GroundStateEnergyCalculator {
	virtual double computeGroundStateEnergy(CompositeSpinInstruction& inst,
			const int nQubits);

	virtual const std::string name() const {
		return "vqe-arpack";
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
#endif /* TASK_TASKS_ARPACK_ARPACKPPGROUNDSTATEENERGYCALCULATOR_HPP_ */
