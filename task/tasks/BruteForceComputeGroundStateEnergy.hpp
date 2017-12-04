#ifndef VQETASKS_BRUTEFORCECOMPUTEGROUNDSTATEENERGY_HPP_
#define VQETASKS_BRUTEFORCECOMPUTEGROUNDSTATEENERGY_HPP_

#include "VQETask.hpp"
#include "CompositeSpinInstruction.hpp"

namespace xacc {
namespace vqe {

class BruteForceComputeGroundStateEnergy: public VQETask {

public:

	BruteForceComputeGroundStateEnergy() {}

	BruteForceComputeGroundStateEnergy(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return the name of this instance.
	 *
	 * @return name The string name
	 */
	virtual const std::string name() const {
		return "vqe-bf-gse";
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
				"Compute Ground State Energy Options");
		desc->add_options()("vqe-ground-state-calculator", value<std::string>(), "");
		return desc;
	}

	virtual bool handleOptions(variables_map& map) {
		return false;
	}

};

class GroundStateEnergyCalculator : public Identifiable {

public:

	virtual double computeGroundStateEnergy(CompositeSpinInstruction& inst,
			const int nQubits) = 0;
	virtual ~GroundStateEnergyCalculator() {}
};

class EigenMatrixXcdGroundStateCalculator: public GroundStateEnergyCalculator {
public:

	virtual double computeGroundStateEnergy(CompositeSpinInstruction& inst,
			const int nQubits);

	virtual const std::string name() const {
		return "vqe-eigen-gs-calculator";
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
