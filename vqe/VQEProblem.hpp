#ifndef VQE_VQEPROBLEM_HPP_
#define VQE_VQEPROBLEM_HPP_

#include "problem.h"
#include "XACC.hpp"
#include "VQEGateFunction.hpp"

namespace xacc {

namespace vqe {

template<typename T>
class VQEProblem : public cppoptlib::Problem<T> {

protected:

	std::shared_ptr<AcceleratorBuffer> buffer;

	std::vector<Kernel<>> kernels;

	int nParameters;

public:

	using typename cppoptlib::Problem<T>::TVector;

	VQEProblem(std::istream& moleculeKernel) :nParameters(0) {

		auto serviceRegistry = xacc::ServiceRegistry::instance();
		auto runtimeOptions = RuntimeOptions::instance();
		(*runtimeOptions)["compiler"] = "fermion";
		runtimeOptions->insert(std::make_pair("state-preparation", "uccsd"));
		runtimeOptions->insert(std::make_pair("n-electrons", "2"));

		// Create the Accelerator
		auto qpu = xacc::getAccelerator("tnqvm");

		// Create the Program
		xacc::Program program(qpu, moleculeKernel);

		// Start compilation
		program.build();

		// Create a buffer of qubits
		std::string nQbtsStr = (*runtimeOptions)["n-qubits"];
		buffer = qpu->createBuffer("qreg", std::stoi(nQbtsStr));

		// Get the Kernels that were created
		kernels = program.getRuntimeKernels();

		// Set the number of VQE parameters
		nParameters = kernels[0].getNumberOfKernelParameters();
	}

	typename cppoptlib::Problem<T>::TVector initializeParameters() {
		return Eigen::VectorXd::Random(nParameters);
	}


	T value(const TVector& x) {

		std::vector<InstructionParameter> parameters;
		for (int i = 0; i < x.rows(); i++) {
			InstructionParameter p(x(i));
			parameters.push_back(p);
		}

		double sum = 0.0;
		for (int i = 0; i < kernels.size(); i++) {
			std::cout << "EXECUTING " << i << " KERNEL\n";

			double expectationValue = 0.0;
			auto kernel = kernels[i];
			auto vqeFunction = std::dynamic_pointer_cast<VQEGateFunction>(
							kernels[i].getIRFunction());

			if (vqeFunction->nInstructions() > 0) {
				std::cout << "QASM:\n" << vqeFunction->toString("qreg") << "\n";
				kernels[i](buffer, parameters);
				// Get Expectation value
				expectationValue = buffer->getExpectationValueZ();
			} else {
				expectationValue = 1.0;
			}

			std::cout << "Found Expectation - " << expectationValue << "\n";

			sum += vqeFunction->coefficient * expectationValue;

			buffer->resetBuffer();
		}

		return sum;
	}
};
}

}

#endif /* VQE_VQEPROBLEM_HPP_ */
