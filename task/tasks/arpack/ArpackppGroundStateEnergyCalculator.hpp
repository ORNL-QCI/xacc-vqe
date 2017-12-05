#ifndef TASK_TASKS_ARPACK_ARPACKPPGROUNDSTATEENERGYCALCULATOR_HPP_
#define TASK_TASKS_ARPACK_ARPACKPPGROUNDSTATEENERGYCALCULATOR_HPP_

#include "BruteForceComputeGroundStateEnergy.hpp"

namespace xacc {
namespace vqe {

template<class T>
class MyMatrix {
protected:

	Eigen::SparseMatrix<double>& m;

public:

	MyMatrix(Eigen::SparseMatrix<double>& mat) : m(mat) {}

	int ncols() { return 16; }

	void MultMv(T* v, T* w) {
		// Execute m * v

		for (int k = 0; k < m.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it;
					++it) {
				if (it.value() != std::complex<double>(0.0, 0.0)) {
					auto r = it.row();
					auto c = it.col();
					auto val = it.value();
//					std::cout << "Setting: w[" << r << "] = " << val*v[c] << ", (" << val << "*" << v[c] << " for c " << c << ") \n";
					w[r] += val * v[c];
				}
			}
		}

	}
};
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
