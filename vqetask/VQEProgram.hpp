#ifndef PROGRAM_VQEPROGRAM_HPP_
#define PROGRAM_VQEPROGRAM_HPP_

#include "Program.hpp"
#include "XACC.hpp"
#include "Function.hpp"
#include "IRGenerator.hpp"

#include <boost/mpi.hpp>
#include "CountGatesOfTypeVisitor.hpp"
#include "Measure.hpp"

namespace xacc {
namespace vqe {

class VQEProgram: public xacc::Program {

public:

	VQEProgram(std::shared_ptr<Accelerator> acc, const std::string& kernelSrc,
			boost::mpi::communicator& c) :
			Program(acc, kernelSrc), nParameters(0), comm(c) {

	}

	VQEProgram(std::shared_ptr<Accelerator> acc,
			const std::string& kernelSource, const std::string& statePrepSrc,
			boost::mpi::communicator& c) :
			Program(acc, kernelSource), statePrepSource(statePrepSrc), nParameters(
					0), comm(c) {

	}

	boost::mpi::communicator& getCommunicator() {
		return comm;
	}

	virtual void build() {
		bool userProvidedKernels = false;

			// This class only takes kernels
			// represented as Fermion Kernels.
			xacc::setCompiler("fermion");

			auto nKernels = 0;
			size_t nPos = src.find("__qpu__", 0);
			while (nPos != std::string::npos) {
				nKernels++;
				nPos = src.find("__qpu__", nPos + 1);
			}

			// Create a buffer of qubits
			nQubits = std::stoi(xacc::getOption("n-qubits"));

			std::vector<double> coeffs;
			// If nKernels > 1, we have non-fermioncompiler kernels
			// so lets check to see if they provided any coefficients
			if (nKernels > 1) { // && boost::contains(src, "coefficients")) {
				// FIXME ADD HOOK TO SET COMPILER FOR THESE
				xacc::setCompiler("scaffold");
				if (xacc::optionExists("vqe-kernels-compiler")) {
					xacc::setCompiler(xacc::getOption("vqe-kernels-compiler"));

				}
				userProvidedKernels = true;
				accelerator->createBuffer("qreg", nQubits);
			}

			addPreprocessor("fcidump-preprocessor");

			// Start compilation
			Program::build();

			// Get the Kernels that were created
			kernels = getRuntimeKernels();

			if (userProvidedKernels) {
				if (boost::contains(src, "pragma")
						&& boost::contains(src, "coefficient")) {
					std::vector<std::string> lines;
					boost::split(lines, src, boost::is_any_of("\n"));
					int counter = 0;
					for (int i = 0; i < lines.size(); ++i) {
						auto line = lines[i];
						if (boost::contains(line, "#pragma vqe-coefficient")) {
							std::vector<std::string> splitspaces;
							boost::split(splitspaces, line, boost::is_any_of(" "));
							boost::trim(splitspaces[2]);
							coeffs.push_back(std::stod(splitspaces[2]));
							InstructionParameter p(
									std::complex<double>(std::stod(splitspaces[2]),
											0.0));
							kernels[counter].getIRFunction()->addParameter(p);
							counter++;
						}
					}
				}
			}

			statePrep = createStatePreparationCircuit();

			// Set the number of VQE parameters
			nParameters = statePrep->nParameters();

	}

	KernelList<> getVQEKernels() {
		return kernels;
	}

	std::shared_ptr<Function> getStatePreparationCircuit() {
		return statePrep;
	}

	const int getNParameters() {
		return nParameters;
	}

	const int getNQubits() {
		return nQubits;
	}

	const std::string getStatePrepType() {
		return statePrepType;
	}

	std::shared_ptr<Accelerator> getAccelerator() {
		return accelerator;
	}

	virtual ~VQEProgram() {}

protected:

	int nQubits = 0;

	std::string statePrepType = "uccsd";

	std::string statePrepSource = "";

	boost::mpi::communicator comm;

	/**
	 * Reference to the state preparation circuit
	 * represented as XACC IR.
	 */
	std::shared_ptr<Function> statePrep;

	/**
	 * Reference to the compiled XACC
	 * Kernels. These kernels each represent
	 * a term in the spin Hamiltonian (the
	 * Hamiltonian produced from a Jordan-Wigner
	 * or Bravi-Kitaev transformation). Each kernel
	 * amounts to a state preparation circuit
	 * followed by appropriately constructed qubit measurements.
	 */
	KernelList<> kernels;

	/**
	 * The number of parameters in the
	 * state preparation circuit.
	 */
	int nParameters;

	std::shared_ptr<Function> createStatePreparationCircuit() {

		if (!statePrepSource.empty()) {
			xacc::setCompiler("scaffold");
			if (xacc::optionExists("vqe-state-prep-kernel-compiler")) {
				xacc::setCompiler(
						xacc::getOption("vqe-state-prep-kernel-compiler"));
			}

			statePrepType = "custom";

			XACCInfo("GENERATING STATE PREP CIRCUIT\n"+statePrepSource+"\n");
			Program p(accelerator, statePrepSource);
			p.build();

			auto kernel = p.getRuntimeKernels()[0];

			return kernel.getIRFunction();
		} else if (xacc::optionExists("vqe-state-prep-kernel")) {
			auto filename = xacc::getOption("vqe-state-prep-kernel");
			std::ifstream filess(filename);

			xacc::setCompiler("scaffold");
			if (xacc::optionExists("vqe-state-prep-kernel-compiler")) {
				xacc::setCompiler(
						xacc::getOption("vqe-state-prep-kernel-compiler"));
			}

			statePrepType = "custom";

			Program p(accelerator, filess);
			p.build();

			auto kernel = p.getRuntimeKernels()[0];

			return kernel.getIRFunction();
		} else {
			if (xacc::optionExists("state-preparation")) {
				statePrepType = xacc::getOption("state-preparation");
			}

			auto statePrepGenerator = ServiceRegistry::instance()->getService<
					IRGenerator>(statePrepType);
			return statePrepGenerator->generate(
					std::make_shared<AcceleratorBuffer>("", nQubits));
		}
	}

};
}
}

#endif
