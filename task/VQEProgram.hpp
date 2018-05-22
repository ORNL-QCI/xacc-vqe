#ifndef PROGRAM_VQEPROGRAM_HPP_
#define PROGRAM_VQEPROGRAM_HPP_

#include "Program.hpp"
#include "XACC.hpp"
#include "Function.hpp"
#include "IRGenerator.hpp"
#include "PauliOperator.hpp"
#include "FermionToSpinTransformation.hpp"

#include "MPIProvider.hpp"
#include "CountGatesOfTypeVisitor.hpp"
#include "Measure.hpp"

#include "IRProvider.hpp"

#include "unsupported/Eigen/CXX11/Tensor"

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
		xacc::error("Error - you have tried to execute the VQEDummyAccelerator. "
				"Please use a real Accelerator.");
	}
	virtual std::vector<std::shared_ptr<AcceleratorBuffer>> execute(
			std::shared_ptr<AcceleratorBuffer> buffer,
			const std::vector<std::shared_ptr<Function>> functions) {
		xacc::error("Error - you have tried to execute the VQEDummyAccelerator. "
						"Please use a real Accelerator.");
		return std::vector<std::shared_ptr<AcceleratorBuffer>>{};
	}
	virtual std::shared_ptr<AcceleratorBuffer> createBuffer(
				const std::string& varId) {
		xacc::error("Error - you have tried to create an AcceleratorBuffer "
				"with the VQEDummyAccelerator. "
						"Please use a real Accelerator.");
		return std::make_shared<AcceleratorBuffer>("",1);
	}
	virtual std::shared_ptr<AcceleratorBuffer> createBuffer(
			const std::string& varId, const int size) {
		xacc::error("Error - you have tried to create an AcceleratorBuffer "
				"with the VQEDummyAccelerator. "
						"Please use a real Accelerator.");
		return std::make_shared<AcceleratorBuffer>("",1);
	}
	virtual bool isValidBufferSize(const int NBits) {
		return false;
	}
	virtual const std::string name() const { return "vqe-dummy"; }
	virtual const std::string description() const {return "";}
};

class VQEProgram: public xacc::Program, public OptionsProvider {

public:

	VQEProgram(std::shared_ptr<Accelerator> acc, PauliOperator& op,
			std::shared_ptr<xacc::Function> sprep,
			std::shared_ptr<Communicator> c) :
			Program(acc, ""), pauli(op), comm(c), nParameters(sprep ? sprep->nParameters() : 0), statePrep(
					sprep), kernels(acc) {
	}

	VQEProgram(std::shared_ptr<Accelerator> acc, const std::string& kernelSource,
			std::shared_ptr<xacc::Function> sprep,
			std::shared_ptr<Communicator> c) :
			Program(acc, kernelSource), comm(c), nParameters(sprep ? sprep->nParameters() : 0), statePrep(
					sprep), kernels(acc) {
	}

	VQEProgram(std::shared_ptr<Accelerator> acc, const std::string& kernelSrc,
			std::shared_ptr<Communicator> c) :
			Program(acc, kernelSrc), nParameters(0), comm(c) {
	}

	VQEProgram(std::shared_ptr<Accelerator> acc,
			const std::string& kernelSource, const std::string& statePrepSrc,
			std::shared_ptr<Communicator> c) :
			Program(acc, kernelSource), statePrepSource(statePrepSrc), nParameters(
					0), comm(c) {
	}

	/**
	 * Return a Boost options_description instance that
	 * describes the options available for this
	 * derived subclass.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		auto desc = std::make_shared<options_description>(
				"VQE Program Options");
		desc->add_options()("correct-readout-errors", "Turn on readout-error correction.");
		return desc;

	}

	/**
	 * Given user-input command line options, perform
	 * some operation. Returns true if runtime should exit,
	 * false otherwise.
	 *
	 * @param map The mapping of options to values
	 * @return exit True if exit, false otherwise
	 */
	virtual bool handleOptions(variables_map& map) {
		return false;
	}

	std::shared_ptr<Communicator> getCommunicator() {
		return comm;
	}

	virtual void build() {

		if (pauli == PauliOperator()) {
			bool userProvidedKernels = false;

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
				if (xacc::optionExists("compiler")) {
					xacc::info("Overridding default compiler to " + xacc::getOption("compiler"));
					xacc::setCompiler(xacc::getOption("compiler"));
				} else {
					xacc::setCompiler("scaffold");
				}
				if (!xacc::optionExists("n-qubits")) {
					xacc::error("You must provide --n-qubits arg if "
							"running with custom hamiltonian kernels.");
				}
				nQubits = std::stoi(xacc::getOption("n-qubits"));
				userProvidedKernels = true;
				accelerator->createBuffer("qreg", nQubits);
			} else {
				xacc::setCompiler("fermion");
			}

			addPreprocessor("fcidump-preprocessor");

			if (xacc::optionExists("correct-readout-errors")) {
				addIRPreprocessor("readout-error-preprocessor");
			}

			// Start compilation
			Program::build();

			if (!userProvidedKernels) {
				std::shared_ptr<FermionToSpinTransformation> transform;
				if (xacc::optionExists("fermion-transformation")) {
					auto transformStr = xacc::getOption(
							"fermion-transformation");
					transform = xacc::getService<
							FermionToSpinTransformation>(transformStr);
				} else {
					transform = xacc::getService<
							FermionToSpinTransformation>("jw");
				}
				pauli = transform->getResult();

				// Rerun the build and get reference to the
				// generated fermionkernel
				xacc::setOption("no-fermion-transformation","");
				auto c = getCompiler("fermion");
				auto ir = c->compile(src, accelerator);
				fermionKernel = std::dynamic_pointer_cast<FermionKernel>(ir->getKernels()[0]);
				xacc::unsetOption("no-fermion-transformation");
			}

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

					pauli.fromXACCIR(xaccIR);
				}
			}

		} else {

			if (xacc::optionExists("correct-readout-errors")) {
				xacc::info("Adding Readout Error IR Preprocessor");
				addIRPreprocessor("readout-error-preprocessor");
			}

			nQubits = std::stoi(xacc::getOption("n-qubits"));
			auto tmpKernels = pauli.toXACCIR()->getKernels();
			xaccIR = xacc::getService<IRProvider>("gate")->createIR();
			for (auto t : tmpKernels) {
				xaccIR->addKernel(t);
			}

			// Execute hardware dependent IR Transformations
			auto accTransforms = accelerator->getIRTransformations();
			for (auto t : accTransforms) {
				xaccIR = t->transform(xaccIR);
			}

			for (auto irp : irpreprocessors) {
				bufferPostprocessors.push_back(irp->process(*xaccIR));
			}

			kernels = getRuntimeKernels();
		}

		// We don't need state prep if we are diagonalizing, profiling,
		// generating openfermion scripts, or if we've already been given
		// an ansatz
		auto task = xacc::getOption("vqe-task");
		if (!statePrep && (task == "vqe" || task == "compute-energy")) {
			xacc::info("Creating a StatePreparation Circuit");
			statePrep = createStatePreparationCircuit();

			// Set the number of VQE parameters
			nParameters = statePrep->nParameters();
		}

		if (statePrep) {

			if (xacc::optionExists("qubit-map")) {
				for (auto pp : irpreprocessors) {
					if (pp->name() == "qubit-map-preprocessor") {

						// Make sure we haven't already mapped it
						auto mapStr = xacc::getOption("qubit-map");
						std::vector<std::string> split;
						boost::split(split, mapStr, boost::is_any_of(","));
						std::vector<int> qubitMap;
						int counter = 0;
						for (auto s : split) {
							auto idx = std::stoi(s);
							qubitMap.push_back(idx);
						}

						bool mapStatePrep = true;
						for (auto& inst : statePrep->getInstructions()) {
							for (auto bit : inst->bits()) {
								if (std::find(qubitMap.begin(), qubitMap.end(),
										bit) != qubitMap.end()) {
									mapStatePrep = false;
									break;
								}
							}
						}

						if (mapStatePrep) {
							auto ir = xacc::getService<IRProvider>("gate")->createIR();
							ir->addKernel(statePrep);
							pp->process(*ir);
						}
					}
				}
			}
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

	void setStatePreparationCircuit(std::shared_ptr<Function> s) {
		statePrep = s;
	}
	const int getNParameters() {
		return nParameters;
	}

	const int getNQubits() {
		return nQubits;
	}

	void setNQubits(const int n) {nQubits = n;}

	const std::string getStatePrepType() {
		return statePrepType;
	}

	std::shared_ptr<Accelerator> getAccelerator() {
		return accelerator;
	}

	const double h_nuclear() {
		if (!fermionKernel) {
			xacc::error("Cannot get E_nuc if you did not compile with FermionCompiler");
		}
		return fermionKernel->E_nuc();
	}

	Eigen::Tensor<std::complex<double>, 2> hpq() {
		if (!fermionKernel) {
			xacc::error("Cannot get h_pq if you did not compile with FermionCompiler");
		}
		return fermionKernel->hpq(nQubits);
	}

	Eigen::Tensor<std::complex<double>, 4> hpqrs() {
		if (!fermionKernel) {
			xacc::error("Cannot get h_pqrs if you did not compile with FermionCompiler");
		}
		return fermionKernel->hpqrs(
				nQubits);
	}

	virtual ~VQEProgram() {
	}

protected:

	int nQubits = 0;

	std::string statePrepType = "uccsd";

	std::string statePrepSource = "";

	std::shared_ptr<Communicator> comm;

	std::shared_ptr<FermionKernel> fermionKernel;

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
			if (xacc::optionExists("compiler")) {
				xacc::setCompiler(
						xacc::getOption("compiler"));
			} else {
				xacc::setCompiler("scaffold");
			}

			statePrepType = "custom";

			Program p(accelerator, statePrepSource);
			p.build();

			auto kernel = p.getRuntimeKernels()[0];

			return kernel.getIRFunction();
		} else if (xacc::optionExists("vqe-ansatz")) {
			auto filename = xacc::getOption("vqe-ansatz");
			std::ifstream filess(filename);

			if (xacc::optionExists("compiler")) {
				xacc::setCompiler(
						xacc::getOption("compiler"));
			} else {
				xacc::setCompiler("scaffold");
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

			auto statePrepGenerator = xacc::getService<
					IRGenerator>(statePrepType);
			return statePrepGenerator->generate(
					std::make_shared<AcceleratorBuffer>("", nQubits));
		}
	}

};
}
}

#endif
