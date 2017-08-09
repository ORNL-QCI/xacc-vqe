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
		runtimeOptions->insert(std::make_pair("n-electrons", "2"));

		auto qpu = xacc::getAccelerator("simple");

		xacc::Program program(qpu, moleculeKernel);

		program.build();

		auto fermionCompiler = serviceRegistry->getService<xacc::Compiler>(
				"fermion");
		int nQubits = fermionCompiler->getNQubitsUsed();

		buffer = qpu->createBuffer("qreg", nQubits);


		runtimeOptions->insert(std::make_pair("n-qubits", std::to_string(nQubits)));
		auto statePrep = serviceRegistry->getService<xacc::IRTransformation>(
				"uccsd");

		program.transformIR(statePrep);

		kernels = program.getRuntimeKernels();
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
			kernels[i](buffer, parameters);

			// Get Expectation value
			double expectationValue = buffer->getAverage();// get from buffer

			auto vqeFunction = std::dynamic_pointer_cast<VQEGateFunction>(
					kernels[i].getIRFunction());
			sum += vqeFunction->coefficient * expectationValue;
		}

		return sum;
	}
};
}

}

#endif /* VQE_VQEPROBLEM_HPP_ */
