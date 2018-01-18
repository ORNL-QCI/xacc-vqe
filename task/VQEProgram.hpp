#ifndef PROGRAM_VQEPROGRAM_HPP_
#define PROGRAM_VQEPROGRAM_HPP_

#include "Program.hpp"
#include "XACC.hpp"
#include "Function.hpp"
#include "IRGenerator.hpp"
#include "PauliOperator.hpp"
#include "FermionToSpinTransformation.hpp"

#include <boost/mpi.hpp>
#include "CountGatesOfTypeVisitor.hpp"
#include "Measure.hpp"

namespace xacc {
namespace vqe {

/**
 * When users just want to run the brute force energy calculator
 * or generate Hamiltonian profiles, they don't actually need
 * a real Accelerator. So here we provide a dummy one just in case,
 * this simplifies the build process for xacc-vqe.
 */
class VQEDummyAccelerator : public xacc::Accelerator {
public:
	virtual void initialize() {}
	virtual AcceleratorType getType() { return AcceleratorType::qpu_gate; }
	virtual std::vector<std::shared_ptr<IRTransformation>> getIRTransformations() {
		return std::vector<std::shared_ptr<IRTransformation>> {};
	}
	virtual void execute(std::shared_ptr<AcceleratorBuffer> buffer,
				const std::shared_ptr<Function> function) {
		XACCError("Error - you have tried to execute the VQEDummyAccelerator. "
				"Please use a real Accelerator.");
	}
	virtual std::vector<std::shared_ptr<AcceleratorBuffer>> execute(
			std::shared_ptr<AcceleratorBuffer> buffer,
			const std::vector<std::shared_ptr<Function>> functions) {
		XACCError("Error - you have tried to execute the VQEDummyAccelerator. "
						"Please use a real Accelerator.");
	}
	virtual std::shared_ptr<AcceleratorBuffer> createBuffer(
				const std::string& varId) {
		XACCError("Error - you have tried to create an AcceleratorBuffer "
				"with the VQEDummyAccelerator. "
						"Please use a real Accelerator.");
	}
	virtual std::shared_ptr<AcceleratorBuffer> createBuffer(
			const std::string& varId, const int size) {
		XACCError("Error - you have tried to create an AcceleratorBuffer "
				"with the VQEDummyAccelerator. "
						"Please use a real Accelerator.");
	}
	virtual bool isValidBufferSize(const int NBits) {
		return false;
	}
	virtual const std::string name() const { return "vqe-dummy"; }
	virtual const std::string description() const {return "";}
};

class VQEProgram: public xacc::Program {

public:

	VQEProgram(std::shared_ptr<Accelerator> acc, PauliOperator& op,
			std::shared_ptr<xacc::quantum::GateFunction> sprep,
			boost::mpi::communicator& c) :
			Program(acc, ""), pauli(op), comm(c), nParameters(sprep ? sprep->nParameters() : 0), statePrep(
					sprep), kernels(acc) {
	}

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

		if (pauli == PauliOperator()) {
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

			std::vector<double> coeffs;
			// If nKernels > 1, we have non-fermioncompiler kernels
			// so lets check to see if they provided any coefficients
			if (nKernels > 1) { // && boost::contains(src, "coefficients")) {
				xacc::setCompiler("scaffold");
				if (xacc::optionExists("vqe-kernels-compiler")) {
					xacc::setCompiler(xacc::getOption("vqe-kernels-compiler"));
				}
				if (!xacc::optionExists("n-qubits")) {
					XACCError("You must provide --n-qubits arg if "
							"running with custom hamiltonian kernels.");
				}
				nQubits = std::stoi(xacc::getOption("n-qubits"));
				userProvidedKernels = true;
				accelerator->createBuffer("qreg", nQubits);
			}

			addPreprocessor("fcidump-preprocessor");

			// Start compilation
			Program::build();

			std::shared_ptr<FermionToSpinTransformation> transform;
			if (xacc::optionExists("fermion-transformation")) {
				auto transformStr = xacc::getOption("fermion-transformation");
				transform = ServiceRegistry::instance()->getService<FermionToSpinTransformation>(
						transformStr);
			} else {
				transform = ServiceRegistry::instance()->getService<FermionToSpinTransformation>(
						"jw");
			}

			pauli = transform->getResult();

			nQubits = std::stoi(xacc::getOption("n-qubits"));

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
							boost::split(splitspaces, line,
									boost::is_any_of(" "));
							boost::trim(splitspaces[2]);
							coeffs.push_back(std::stod(splitspaces[2]));
							InstructionParameter p(
									std::complex<double>(
											std::stod(splitspaces[2]), 0.0));
							InstructionParameter q(
									(kernels[counter].getIRFunction()->nInstructions()
											== 0 ? 1 : 0));
							kernels[counter].getIRFunction()->addParameter(p);
							kernels[counter].getIRFunction()->addParameter(q);
							counter++;
						}
					}
				}
			}

		} else {
			nQubits = std::stoi(xacc::getOption("n-qubits"));
			auto tmpKernels = pauli.toXACCIR()->getKernels();
			for (auto t : tmpKernels) {
				kernels.push_back(Kernel<>(accelerator,t));
			}
		}

		// We don't need state prep if we are brute
		// force computing ground state energy
		if (!statePrep && xacc::getOption("vqe-task") != "vqe-diagonalize"
				&& xacc::getOption("vqe-task") != "vqe-profile"
				&& xacc::getOption("vqe-task")
						!= "vqe-openfermion-eigenspectrum") {
			XACCInfo("Creating a StatePreparation Circuit");
			statePrep = createStatePreparationCircuit();

			// Set the number of VQE parameters
			nParameters = statePrep->nParameters();
		}

	}

	PauliOperator getPauliOperator() {
		return pauli;
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

	virtual ~VQEProgram() {
	}

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

	PauliOperator pauli = PauliOperator();

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
