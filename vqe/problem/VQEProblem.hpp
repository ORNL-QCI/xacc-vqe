#ifndef VQE_VQEPROBLEM_HPP_
#define VQE_VQEPROBLEM_HPP_

#include "problem.h"
#include "XACC.hpp"
#include <boost/math/constants/constants.hpp>

#include "GateQIR.hpp"

#include "InstructionIterator.hpp"

#include <boost/python.hpp>
using namespace boost::python;
using namespace xacc::quantum;

namespace xacc {

namespace vqe {

class VQEProblem : public cppoptlib::Problem<double> {

protected:

	std::shared_ptr<Accelerator> qpu;

	std::shared_ptr<Program> program;

	std::vector<Kernel<>> kernels;

	int nParameters;

	int nQubits;

	std::shared_ptr<Function> statePrep;

public:

	double currentEnergy;

	VQEProblem(std::istream& moleculeKernel) : nParameters(0), currentEnergy(0.0) {
		xacc::setCompiler("fermion");
		xacc::setOption("n-electrons", "2");

		// Create the Accelerator
		auto qpu = xacc::getAccelerator("tnqvm");

		// Create the Program
		program = std::make_shared<Program>(qpu, moleculeKernel);

		// Start compilation
		program->build();

		// Create a buffer of qubits
		nQubits = std::stoi(xacc::getOption("n-qubits"));

		// Get the Kernels that were created
		kernels = program->getRuntimeKernels();

		statePrep = createStatePreparationCircuit();

		// Set the number of VQE parameters
		nParameters = statePrep->nParameters();

		if (xacc::optionExists("vqe-print-scaffold-source")) {
			printScaffoldSourceCode();
		}
	}

	Eigen::VectorXd initializeParameters() {
		std::srand(time(0));
		auto pi = boost::math::constants::pi<double>();
		// Random parameters between -pi and pi
		auto rand = -1.0 * pi * Eigen::VectorXd::Ones(nParameters)
				+ (Eigen::VectorXd::Random(nParameters) * 0.5
						+ Eigen::VectorXd::Ones(nParameters) * 0.5)
						* (pi - (-1 * pi));
		return rand;
	}


	double value(const Eigen::VectorXd& x) {

		auto evaluatedStatePrep = evaluateStatePreparationCircuit(x);

		double sum = 0.0, localExpectationValue = 0.0;
#pragma omp parallel for reduction (+:sum)
		for (int i = 0; i < kernels.size(); i++) {

			// Get the ith Kernel
			auto irFunction = kernels[i].getIRFunction();
			irFunction->insertInstruction(0, evaluatedStatePrep);
			double coeff = std::real(
					boost::get<std::complex<double>>(
							irFunction->getParameter(0)));

			auto localqpu = xacc::getAccelerator("tnqvm");
			auto buff = localqpu->createBuffer("qreg", nQubits);
			localqpu->execute(buff, irFunction);

			// Get Expectation value
			if (boost::get<int>(irFunction->getParameter(1))) {
				localExpectationValue = 1.0;
			} else {
				localExpectationValue = buff->getExpectationValueZ();
			}

			// Sum up the expectation values
			sum += coeff * localExpectationValue;
			irFunction->removeInstruction(0);
		}

		currentEnergy = sum;

		std::stringstream ss;
		ss << x.transpose();
		XACCInfo("Computed VQE Energy = " + std::to_string(sum) + " at (" + ss.str() + ")");
		return sum;
	}

	class VQECriteria : public cppoptlib::Criteria<double> {
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

private:

	std::shared_ptr<Function> createStatePreparationCircuit() {
		// Create the State Preparation circuit
		auto tempIR = std::make_shared<GateQIR>();
		tempIR->addKernel(
				std::make_shared<GateFunction>("temp",
						std::vector<InstructionParameter> {
								InstructionParameter(std::complex<double>(1.0)),
								InstructionParameter(0) }));
		auto statePrepIRTransform = ServiceRegistry::instance()->getService<
				IRTransformation>(
				xacc::optionExists("state-preparation") ?
						xacc::getOption("state-preparation") : "uccsd");
		auto statePrepIR = statePrepIRTransform->transform(tempIR);
		return std::dynamic_pointer_cast<Function>(
				statePrepIR->getKernels()[0]->getInstruction(0));
	}

	std::shared_ptr<Function> evaluateStatePreparationCircuit(const TVector& x) {

		auto pi = boost::math::constants::pi<double>();
		std::vector<std::string> variableNames;
		for (int i = 0; i < x.rows(); i++) {
			variableNames.push_back(
					boost::get<std::string>(statePrep->getParameter(i)));
		}

		Py_Initialize();
		const std::string code = R"code(val = eval(exprStr))code";
		// Retrieve the main module.
		object main = import("__main__");
		object main_namespace = main.attr("__dict__");

		auto evaluatedStatePrep = std::make_shared<GateFunction>("stateprep");
		for (auto inst : statePrep->getInstructions()) {
			if (!inst->isComposite() && inst->isParameterized()
					&& inst->getParameter(0).which() == 3) {
				int idx = -1;
				auto expression = boost::get<std::string>(inst->getParameter(0));
				for (int i = 0; i < nParameters; i++) {
					if (boost::contains(expression, variableNames[i])) {
						idx = i;
					}
				}
				std::stringstream ss1;
				ss1 << x(idx);
				boost::replace_all(expression,
						variableNames[idx], ss1.str());
				main_namespace["exprStr"] = expression;
				exec(code.c_str(), main_namespace);
				double res = extract<double>(main_namespace["val"]);
				if (res < 0.0) {
					res = 4 * pi + res;
				}
				InstructionParameter p(res);
				auto updatedInst = GateInstructionRegistry::instance()->create(
						inst->getName(), inst->bits());
				updatedInst->setParameter(0, p);
				evaluatedStatePrep->addInstruction(updatedInst);
			} else {
				evaluatedStatePrep->addInstruction(inst);
			}
		}

		Py_Finalize();

		return evaluatedStatePrep;

	}

	void printScaffoldSourceCode() {
		int counter = 0;
		for (auto k : kernels) {
			auto scaffold = xacc::getCompiler("scaffold");
			auto f = k.getIRFunction();
			f->insertInstruction(0, statePrep);
			auto srcStr = scaffold->translate("qreg", f);
			f->removeInstruction(0);

			boost::filesystem::path dir("scaffold_source");
			if (!boost::filesystem::exists(dir)) {
				if (!boost::filesystem::create_directory(dir)) {
					XACCError("Could not create scaffold_source directory.");
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
};
}

}

#endif
//
//auto scaffold = xacc::getCompiler("scaffold");
//	auto f = kernel.getIRFunction();
//	auto srcStr = scaffold->translate("qreg", f);
//
//	boost::filesystem::path dir("temp");
//	if (!boost::filesystem::exists(dir)) {
//		if (!boost::filesystem::create_directory(dir)) {
//			XACCError("Could not create scaffold_source directory.");
//		}
//	}
//	std::ofstream out("temp/kernel_" + std::to_string(i) + "_at_"
//					+ std::to_string(boost::get<double>(parameters[0])) + "_" + std::to_string(boost::get<double>(parameters[1]))
//					+ ".hpp");
//	out << srcStr;
//	out.flush();
//	out.close();

/*
			auto statePrep = std::dynamic_pointer_cast<Function>(k.getIRFunction()->getInstruction(0));
			InstructionIterator it2(statePrep);
			it2.next();
			int counter = 0, index = 0;
			while (it2.hasNext()) {
				auto next = it2.next();
				if (next->getName() == "Rz"
						&& std::fabs(boost::get<double>(next->getParameter(0)))
								< 1e-12) {
//					std::cout << "Found an Rz with 0 angle\n";

//					next->disable();
					// If I'm here, I have counter instructions in front of me
					// So disable the counter in front and counter behind
					while(counter != 0) {
//					for (int i = index-counter; i <= 2*(counter+1); i++) {
//						std::cout << counter << ", " << index << " | Disabling " << statePrep->getInstruction(index-counter)->toString("qreg") << " and " << statePrep->getInstruction(index+counter)->toString("qreg") << "\n";
//						statePrep->getInstruction(index-counter)->disable();
//						statePrep->getInstruction(index+counter)->disable();
						counter--;
					}

					// reset the counter
//					counter = 0;
				}

				if (next->getName() != "X" && next->isEnabled()) {
					counter++;
				}

				index++;
			}*/
