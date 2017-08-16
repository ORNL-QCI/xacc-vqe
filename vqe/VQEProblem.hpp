#ifndef VQE_VQEPROBLEM_HPP_
#define VQE_VQEPROBLEM_HPP_

#include "problem.h"
#include "XACC.hpp"
#include <boost/math/constants/constants.hpp>

#include "VQEGateFunction.hpp"
#include "InstructionIterator.hpp"

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

	double currentEnergy;

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

		if (xacc::optionExists("vqe-print-scaffold-source")) {
			int counter = 0;
			for (auto k : kernels) {
				auto scaffold = xacc::getCompiler("scaffold");
				auto f = k.getIRFunction();
				auto srcStr = scaffold->translate("qreg", f);

				boost::filesystem::path dir("scaffold_source");
				if (!boost::filesystem::exists(dir)) {
					if (!boost::filesystem::create_directory(dir)) {
						XACCError(
								"Could not create scaffold_source directory.");
					}
				}
				std::ofstream out(
						"scaffold_source/"
								+ xacc::getOption("vqe-print-scaffold-source")
								+ "_kernel_" + std::to_string(counter) + ".hpp");
				out << srcStr;
				out.flush();
				out.close();
				counter++;
			}

			if (xacc::optionExists("vqe-exit-after-scaffold")) {
				xacc::Finalize();
				exit(0);
			}
		}
	}

	typename cppoptlib::Problem<T>::TVector initializeParameters() {
		std::srand(time(0));
		auto rand = Eigen::VectorXd::Random(nParameters);
		std::cout << "InitialParams: " << rand.transpose() << "\n";
		return rand;
	}


	T value(const TVector& x) {
		auto pi = boost::math::constants::pi<double>();
		std::vector<InstructionParameter> parameters;
		for (int i = 0; i < x.rows(); i++) {
			InstructionParameter p(x(i));
			parameters.push_back(p);
		}

//		InstructionParameter p1(0.0);
//		InstructionParameter p2(.05677);
//		parameters.clear();
//		parameters.push_back(p1);
//		parameters.push_back(p2);

		// Evaluate all parameters first,
		// since this invokes the Python Interpreter (for now)
		for (auto k : kernels) {
			k.evaluateParameters(parameters);
			InstructionIterator it(k.getIRFunction());
			while (it.hasNext()) {
				// Get the next node in the tree
				auto nextInst = it.next();
				auto gateName = nextInst->getName();
				if (gateName == "Rx" || gateName == "Rz") {
					auto angle =
							boost::get<double>(nextInst->getParameter(0));
					if (angle < 0.0) {
						InstructionParameter newParam(
								gateName == "Rx" ?
										(std::fabs(angle) + 3 * pi) :
										(std::fabs(angle) + 4 * pi));
						nextInst->setParameter(0, newParam);
					}
				}
			}

		}

		double sum = 0.0, localExpectationValue = 0.0;
#pragma omp parallel for reduction (+:sum)
		for (int i = 0; i < kernels.size(); i++) {

			// Get the ith Kernel
			auto kernel = kernels[i];

			auto scaffold = xacc::getCompiler("scaffold");
			auto f = kernel.getIRFunction();
			auto srcStr = scaffold->translate("qreg", f);

			boost::filesystem::path dir("temp");
			if (!boost::filesystem::exists(dir)) {
				if (!boost::filesystem::create_directory(dir)) {
					XACCError("Could not create scaffold_source directory.");
				}
			}
			std::ofstream out("temp/kernel_" + std::to_string(i) + "_at_"
							+ std::to_string(boost::get<double>(parameters[0])) + "_" + std::to_string(boost::get<double>(parameters[1]))
							+ ".hpp");
			out << srcStr;
			out.flush();
			out.close();

			// We need the reference to the IR Function
			// in order to get the leading coefficient
			auto vqeFunction = std::dynamic_pointer_cast<VQEGateFunction>(
							kernel.getIRFunction());
			double coeff = vqeFunction->coefficient;

			// We only need a temporary AcceleratorBuffer,
			// so just create it here for this thread's iteration
			auto buff = qpu->createBuffer("qreg", nQubits);

			// If we have instructions, execute the kernel
			// Execute!
			if (vqeFunction->nInstructions() > 2) kernel(buff);

			// Get Expectation value
			localExpectationValue = buff->getExpectationValueZ();

			// Sum up the expectation values
			sum += coeff * localExpectationValue;
			 std::cout << " Kernel " << i << " Expectation = " << localExpectationValue << ", coeff = " << coeff << ": " << sum << "\n";

		}

		currentEnergy = sum;

		std::stringstream ss;
		ss << parameters[0] << " " << parameters[1];//x.transpose();
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
