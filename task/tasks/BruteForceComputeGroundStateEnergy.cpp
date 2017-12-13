#include "BruteForceComputeGroundStateEnergy.hpp"
#include "FermionToSpinTransformation.hpp"
#include <iostream>
#include <unordered_map>
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
//	= inst.toSparseRealMatrix(
//			nQubits);

	auto nTerms = inst.nInstructions();

	// Get Identity coefficient
	std::complex<double> identityCoeff(0.0, 0.0);
	for (int i = 0; i < nTerms; i++) {
		std::shared_ptr<SpinInstruction> spinInst = std::dynamic_pointer_cast<
				SpinInstruction>(inst.getInstruction(i));
		if (spinInst->isIdentity()) {
			identityCoeff = spinInst->coefficient;
			break;
		}
	}

	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	// Generate all bit strings
	std::vector<std::string> bitStrings(dim);
#pragma omp parallel for
	for (std::uint64_t j = 0; j < dim; j++) {

		std::stringstream s;
		for (int k = nQubits - 1; k >= 0; k--) {
			s << ((j >> k) & 1);
		}
		bitStrings[j] = s.str();
	}

	int Istart =0;
	int Iend = dim;

	using Triplet = Eigen::Triplet<double>;
	using IndexPair = std::pair<std::uint64_t, std::uint64_t>;

	std::unordered_map<IndexPair,
					std::complex<double>,
					boost::hash<IndexPair> > nonZeros;

	for (int i = 0; i < dim; i++)
			nonZeros.insert(
					std::make_pair(IndexPair { i, i },
							std::complex<double>(0, 0)));

	std::vector<Triplet> triplets;
	for (std::uint64_t myRow = Istart; myRow < Iend; myRow++) {

		XACCInfo(
				"Matrix Construction for rank " + std::to_string(0)
						+ ", row " + std::to_string(myRow));

		if (identityCoeff != std::complex<double>(0.0, 0.0)) {
			nonZeros[std::make_pair(myRow,myRow)] += identityCoeff;
		}

		for (int i = 0; i < nTerms; i++) {

			std::shared_ptr<SpinInstruction> spinInst =
					std::dynamic_pointer_cast<SpinInstruction>(
							inst.getInstruction(i));
			std::pair<std::string, std::complex<double>> newBitStrCoeff;

			if (!spinInst->isIdentity()) {
				if (spinInst->isDiagonal()) {
					newBitStrCoeff = spinInst->computeActionOnBra(
							bitStrings[myRow]);
					nonZeros[std::make_pair(myRow,myRow)] += newBitStrCoeff.second;
				} else {
					newBitStrCoeff = spinInst->computeActionOnBra(
							bitStrings[myRow]);
					std::uint64_t k = std::stol(newBitStrCoeff.first, nullptr,
							2);

					auto p = std::make_pair(myRow, k);
					if (nonZeros.find(p) != nonZeros.end()) {
						nonZeros[p] += newBitStrCoeff.second;
					} else {
						nonZeros.insert(
								std::make_pair(p, newBitStrCoeff.second));
					}
				}
			}
		}
	}

	for (auto& kv : nonZeros) {
		triplets.push_back(Triplet(kv.first.first, kv.first.second, std::real(kv.second)));
	}

	Eigen::SparseMatrix<double> ham(dim,dim);
	ham.setFromTriplets(triplets.begin(), triplets.end());
	ham.makeCompressed();

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
