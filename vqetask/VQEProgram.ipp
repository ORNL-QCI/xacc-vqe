#include "VQEProgram.hpp"
namespace xacc {
namespace vqe {

void VQEProgram::build() {

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
		accelerator->createBuffer("qreg", 2);
	}

	addPreprocessor("fcidump-preprocessor");

	// Start compilation
	Program::build();

	// Create a buffer of qubits
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

	if (xacc::optionExists("vqe-print-stats")) {

		if (comm.rank() == 0) {
			XACCInfo("Number of Qubits = " + xacc::getOption("n-qubits"));
			XACCInfo(
					"Number of Hamiltonian Terms = "
							+ std::to_string(kernels.size()));
			XACCInfo("State Prep Type: " + statePrepType);
			XACCInfo(
					"Number of Variational Parameters = "
							+ std::to_string(nParameters));
			std::unordered_map<int, int> countKLocals;
			for (auto k : kernels) {
				// There is a measurement for each of the pauli's
				// in a given hamiltonian term, so just count them
				xacc::quantum::CountGatesOfTypeVisitor<xacc::quantum::Measure> visitor(
						k.getIRFunction());
				auto kLocalCount = visitor.countGates();
				auto search = countKLocals.find(kLocalCount);
				if (search != countKLocals.end()) {
					search->second++;
				} else {
					countKLocals[kLocalCount] = 1;
				}
			}

			std::map<int, int> ordered(countKLocals.begin(),
					countKLocals.end());
			for (auto& it : ordered) {
				XACCInfo(
						"N k-Local Terms (k,N) = (" + std::to_string(it.first)
								+ ", " + std::to_string(it.second) + ")");
			}
		}

		exit(0);
	}
}
/*
std::shared_ptr<Function> VQEProgram::createStatePreparationCircuit() {

	if (xacc::optionExists("vqe-state-prep-kernel")) {
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
	} else if (!statePrepSource.empty()) {
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
	} else {
		if (xacc::optionExists("state-preparation")) {
			statePrepType = xacc::getOption("state-preparation");
		}

		auto statePrepGenerator = ServiceRegistry::instance()->getService<
				IRGenerator>(statePrepType);
		return statePrepGenerator->generate(
				std::make_shared<AcceleratorBuffer>("", nQubits));
	}
}*/

}
}
