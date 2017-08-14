#ifndef VQE_VQEPROBLEM_HPP_
#define VQE_VQEPROBLEM_HPP_

#include <omp.h>
#include "problem.h"
#include "XACC.hpp"

#include "VQEGateFunction.hpp"

namespace xacc {

namespace vqe {

template<typename T>
class VQEProblem : public cppoptlib::Problem<T> {

protected:

	std::shared_ptr<Accelerator> qpu;

	std::vector<Kernel<>> kernels;

	int nParameters;

	int nQubits;

public:

	using typename cppoptlib::Problem<T>::TVector;

	VQEProblem(std::istream& moleculeKernel) : nParameters(0), currentEnergy(0.0) {
		xacc::setCompiler("fermion");
		xacc::setOption("state-preparation", "uccsd");
		xacc::setOption("n-electrons", "2");

		// Create the Accelerator
		qpu = xacc::getAccelerator("tnqvm");

		// Create the Program
		xacc::Program program(qpu, moleculeKernel);

		// Start compilation
		program.build();

		// Create a buffer of qubits
		nQubits = std::stoi(xacc::getOption("n-qubits"));

		// Get the Kernels that were created
		kernels = program.getRuntimeKernels();

		// Set the number of VQE parameters
		nParameters = kernels[0].getNumberOfKernelParameters();
	}

	typename cppoptlib::Problem<T>::TVector initializeParameters() {
		std::srand(time(0));
		auto rand = Eigen::VectorXd::Random(nParameters);
		std::cout << "InitialParams: " << rand.transpose() << "\n";
		return rand;
	}

	double currentEnergy;

	T value(const TVector& x) {

		std::vector<InstructionParameter> parameters;
		for (int i = 0; i < x.rows(); i++) {
			InstructionParameter p(x(i));
			parameters.push_back(p);
		}

		// Evaluate all parameters first,
		// since this invokes the Python Interpreter (for now)
		for (auto k : kernels) {
			k.evaluateParameters(parameters);
		}

		double sum = 0.0, localExpectationValue = 0.0;
#pragma omp parallel for reduction (+:sum)
		for (int i = 0; i < kernels.size(); i++) {

			// Get the ith Kernel
			auto kernel = kernels[i];

			// We need the reference to the IR Function
			// in order to get the leading coefficient
			auto vqeFunction = std::dynamic_pointer_cast<VQEGateFunction>(
							kernel.getIRFunction());
			double coeff = vqeFunction->coefficient;

			// We only need a temporary AcceleratorBuffer,
			// so just create it here for this thread's iteration
			auto buff = qpu->createBuffer("qreg", nQubits);

			// If we have instructions, execute the kernel
			// If not just set the local expectation value to 1
			if (vqeFunction->nInstructions() > 0) {

				// Execute!
				kernel(buff);

				// Get Expectation value
				localExpectationValue = buff->getExpectationValueZ();
			} else {
				localExpectationValue = 1.0;
			}

			// Sum up the expectation values
			sum += coeff * localExpectationValue;
		}

		currentEnergy = sum;

		std::stringstream ss;
		ss << x.transpose();
		XACCInfo("Computed VQE Energy = " + std::to_string(sum) + " at (" + ss.str() + ")");
		return sum;
	}

	class VQECriteria : public cppoptlib::Criteria<T> {
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
}

}

#endif


//auto vis = std::make_shared<xacc::vqe::PrintScaffoldVisitor>("qreg");
//
//		auto f = kernels[10].getIRFunction();
//		f->evaluateVariableParameters(parameters);
//
//		InstructionIterator it(f);
//		while (it.hasNext()) {
//			// Get the next node in the tree
//			auto nextInst = it.next();
//			if (nextInst->isEnabled())
//				nextInst->accept(vis);
//		}
//
//		std::cout << "SCAFFOLD SRC:\n" << vis->getScaffoldString() << "\n";
