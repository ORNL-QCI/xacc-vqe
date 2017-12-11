#include "CountGatesOfTypeVisitor.hpp"
#include "Measure.hpp"
#include "VQEProgram.hpp"
#include "FermionToSpinTransformation.hpp"
#include "GenerateOpenFermionEigenspectrumScript.hpp"

namespace xacc {
namespace vqe {


VQETaskResult GenerateOpenFermionEigenspectrumScript::execute(
		Eigen::VectorXd parameters) {

	auto kernels = program->getVQEKernels();
	auto statePrepType = program->getStatePrepType();
	auto comm = program->getCommunicator();
	auto nParameters = program->getNParameters();

	if (comm.rank() == 0) {
		std::string defaultFileName = "gen_openfermion_script.txt";
		if (xacc::optionExists("vqe-openfermion-script-name")) {
			defaultFileName = xacc::getOption("vqe-openfermion-script-name");
		}
		std::ofstream out(defaultFileName);

		std::shared_ptr<IRTransformation> transform;
		if (xacc::optionExists("fermion-transformation")) {
			auto transformStr = xacc::getOption("fermion-transformation");
			transform = ServiceRegistry::instance()->getService<IRTransformation>(
					transformStr);
		} else {
			transform = ServiceRegistry::instance()->getService<IRTransformation>(
					"jordan-wigner");
		}

		auto hamiltonianInstruction = std::dynamic_pointer_cast<
				FermionToSpinTransformation>(transform)->getResult();

		std::stringstream s;

		s << "from openfermion.ops import QubitOperator\n";
		s << "from openfermion.utils import eigenspectrum\n";
		s << "op = \\\n";

		int nInsts = hamiltonianInstruction.nInstructions();
		for (int i = 0; i < nInsts; i++) {
			auto t = hamiltonianInstruction.getInstruction(i);
			auto spin = std::dynamic_pointer_cast<SpinInstruction>(t);
			std::string termStr = "";
			for (auto term : spin->getTerms()) {
				if (term.second != "I") {
					termStr += " " + term.second + std::to_string(term.first) + " ";
				} else {
					termStr += "";
				}
			}

			if (i == nInsts-1) {
				s << "   QubitOperator('" << termStr << "', complex" << spin->coefficient << ")\n";
			} else {
				s << "   QubitOperator('" << termStr << "', complex" << spin->coefficient << ") + \\\n";
			}
		}

		s << "\n\nes = eigenspectrum(op)\n";
		s << "print('Energies = ', es)\n";

		out << s.str();
		out.flush();
		out.close();
	}

	return VQETaskResult{};
}

}
}
