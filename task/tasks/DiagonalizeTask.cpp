#include "FermionToSpinTransformation.hpp"
#include <iostream>
#include <unordered_map>

#include "DiagonalizeTask.hpp"

namespace xacc {
namespace vqe {


VQETaskResult DiagonalizeTask::execute(
		Eigen::VectorXd parameters) {
	int nQubits = std::stoi(xacc::getOption("n-qubits"));

	boost::mpi::communicator world;

	auto hamiltonianInstruction = program->getPauliOperator();

	std::shared_ptr<DiagonalizeBackend> backend;
	if (xacc::optionExists("diagonalize-backend")) {
		auto str = xacc::getOption("diagonalize-backend");
		backend = ServiceRegistry::instance()->getService<
				DiagonalizeBackend>(str);
	} else {
		backend = ServiceRegistry::instance()->getService<
				DiagonalizeBackend>("diagonalize-eigen");
	}

	auto energy = backend->diagonalize(hamiltonianInstruction,
			nQubits);

	VQETaskResult result;
	result.angles = parameters;
	result.energy = energy;

	return result;
}

double EigenDiagonalizeBackend::diagonalize(
		PauliOperator& inst, const int nQubits) {
	boost::mpi::communicator world;
	int rank = world.rank();
	int nRanks = world.size();
	std::complex<double> gsReal;
	auto nTerms = inst.nTerms();

	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	auto getBitStrForIdx = [&](std::uint64_t i) {
		std::stringstream s;
		for (int k = nQubits - 1; k >= 0; k--) s << ((i >> k) & 1);
		return s.str();
	};

	Eigen::MatrixXcd A(dim, dim);
	A.setZero();

	if (rank == 0) xacc::info(
			"Building Matrix for Eigen.");

	for (std::uint64_t myRow = 0; myRow < dim; myRow++) {
		auto rowBitStr = getBitStrForIdx(myRow);
		auto results = inst.computeActionOnBra(rowBitStr);
		for (auto& result : results) {
			std::uint64_t k = std::stol(result.first, nullptr, 2);
			A(myRow, k) += result.second;
		}
	}

	if (rank == 0) xacc::info(
			"Done building Matrix for Eigen.");


	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(
			A);
	auto eigenvalues = es.eigenvalues();
	gsReal = eigenvalues(0);
	std::stringstream ss;
	ss << std::setprecision(12) << gsReal;
	if (rank == 0) xacc::info("Lowest Eigenvalue = " + ss.str());
	return std::real(gsReal);
}

}
}
