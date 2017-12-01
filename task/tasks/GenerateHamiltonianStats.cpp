#include "GenerateHamiltonianStats.hpp"

#include "CountGatesOfTypeVisitor.hpp"
#include "Measure.hpp"
#include "VQEProgram.hpp"

namespace xacc {
namespace vqe {


VQETaskResult GenerateHamiltonianStats::execute(
		Eigen::VectorXd parameters) {

	auto kernels = program->getVQEKernels();
	auto statePrepType = program->getStatePrepType();
	auto comm = program->getCommunicator();
	auto nParameters = program->getNParameters();

	if (comm.rank() == 0) {
		std::ofstream out("gen_hamiltonian_stats_out.txt");

		std::stringstream s;
		s << "Number of Qubits = " << xacc::getOption("n-qubits") << "\n";
		s << "Number of Hamiltonian Terms = "
				<< std::to_string(kernels.size()) << "\n";
		s << "State Prep Type: " << statePrepType << "\n";
		s << "Number of Variational Parameters = "
				<< std::to_string(nParameters) << "\n";

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
			s << "N k-Local Terms (k,N) = (" << std::to_string(it.first)
			<< ", " << std::to_string(it.second) << ")" << "\n";
		}

		out << s.str();
		out.flush();
		out.close();
	}

	return VQETaskResult{};
}

}
}
