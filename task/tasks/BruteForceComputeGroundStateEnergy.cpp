#include "BruteForceComputeGroundStateEnergy.hpp"
#include "FermionToSpinTransformation.hpp"

namespace xacc {
namespace vqe {


VQETaskResult BruteForceComputeGroundStateEnergy::execute(
		Eigen::VectorXd parameters) {
	int nQubits = std::stoi(xacc::getOption("n-qubits"));

	boost::mpi::communicator world;

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

	Eigen::MatrixXcd ham = hamiltonianInstruction.toMatrix(nQubits);

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(ham);
	Eigen::VectorXcd eigenvalues = es.eigenvalues();

	std::stringstream ss;
	ss << eigenvalues.transpose();
	if (world.rank() == 0) XACCInfo("HamiltonianEigenvalues:\n" + ss.str());

	return std::vector<std::pair<Eigen::VectorXd, double>> { { parameters,
			std::real(eigenvalues(0)) } };
}

}
}
