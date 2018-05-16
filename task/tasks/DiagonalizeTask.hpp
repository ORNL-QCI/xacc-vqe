#ifndef VQETASKS_DIAGONALIZETASK_HPP_
#define VQETASKS_DIAGONALIZETASK_HPP_

#include "VQETask.hpp"
#include "PauliOperator.hpp"

namespace xacc {
namespace vqe {

class DiagonalizeTask : public VQETask {

public:

	DiagonalizeTask() {}

	DiagonalizeTask(std::shared_ptr<VQEProgram> prog) :
			VQETask(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters);

	/**
	 * Return the name of this instance.
	 *
	 * @return name The string name
	 */
	virtual const std::string name() const {
		return "vqe-diagonalize";
	}

	/**
	 * Return the description of this instance
	 * @return description The description of this object.
	 */
	virtual const std::string description() const {
		return "Diagonalize the Qubit Hamiltonian "
				"and return the lowest energy eigenvalue.";
	}

	/**
	 * Return an empty options_description, this is for
	 * subclasses to implement.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		auto desc = std::make_shared<options_description>(
				"Diagonalize Options");
		desc->add_options()("diagonalize-backend", value<std::string>(),
							"The backend to use to compute the Hamiltonian eigenspectrum")
			("diag-number-symmetry","Reduce the dimensionality of the problem by considering Hamiltonian subspace spanned by NELEC occupations.");
		return desc;
	}

	virtual bool handleOptions(variables_map& map) {
		return false;
	}

};

class DiagonalizeBackend : public Identifiable {
public:
	virtual double diagonalize(std::shared_ptr<VQEProgram> prog) = 0;
	virtual ~DiagonalizeBackend() {}
};

class EigenDiagonalizeBackend: public DiagonalizeBackend {
public:

	virtual double diagonalize(std::shared_ptr<VQEProgram> prog);

	virtual const std::string name() const {
		return "diagonalize-eigen";
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
