#include "BruteForceComputeGroundStateEnergy.hpp"
#include "FermionToSpinTransformation.hpp"
#include <iostream>

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

	auto calc = ServiceRegistry::instance()->getService<
			GroundStateEnergyCalculator>("vqe-eigen-gs-calculator");

	if (xacc::optionExists("vqe-ground-state-calculator")) {
		auto str = xacc::getOption("vqe-ground-state-calculator");
		calc = ServiceRegistry::instance()->getService<
				GroundStateEnergyCalculator>(str);

	}

	auto energy = calc->computeGroundStateEnergy(hamiltonianInstruction,
			nQubits);
	return std::vector<std::pair<Eigen::VectorXd, double>> { { parameters,
			energy } };
}

double EigenMatrixXcdGroundStateCalculator::computeGroundStateEnergy(
		CompositeSpinInstruction& inst, const int nQubits) {
	boost::mpi::communicator world;
	Eigen::SparseMatrix<double> ham = inst.toSparseRealMatrix(
			nQubits);
	Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(
			ham);
	auto eigenvalues = es.eigenvalues();
	std::stringstream ss;
	ss << eigenvalues.transpose();
	if (world.rank() == 0)
		XACCInfo("HamiltonianEigenvalues:\n" + ss.str());
	return std::real(eigenvalues(0));
}

}
}
