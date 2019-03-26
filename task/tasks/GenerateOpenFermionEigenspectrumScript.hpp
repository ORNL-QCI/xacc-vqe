#ifndef VQETASKS_GENERATEOPENFERMIONSCRIPT_HPP_
#define VQETASKS_GENERATEOPENFERMIONSCRIPT_HPP_

#include "VQETask.hpp"

namespace xacc {
namespace vqe {

class GenerateOpenFermionEigenspectrumScript: public VQETask {

public:

	GenerateOpenFermionEigenspectrumScript() {}

	GenerateOpenFermionEigenspectrumScript(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return the name of this instance.
	 *
	 * @return name The string name
	 */
	virtual const std::string name() const {
		return "vqe-openfermion-eigenspectrum";
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
	virtual OptionPairs getOptions() {
		OptionPairs desc {{"vqe-openfermion-eigenspectrum-script-name", "The name of the file to save."}};
		return desc;
	}

};

}
}
#endif
