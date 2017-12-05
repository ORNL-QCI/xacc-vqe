#include "SlepcGroundStateEnergyCalculator.hpp"

namespace xacc {
namespace vqe {
double SlepcGroundStateEnergyCalculator::computeGroundStateEnergy(
		CompositeSpinInstruction& inst, const int nQubits) {
	boost::mpi::communicator world;

	using SparseComplexMatrix = Eigen::SparseMatrix<std::complex<double>>;
	using SparseRealMatrix = Eigen::SparseMatrix<double>;
	SparseComplexMatrix spMat = inst.toSparseMatrix(nQubits);
	spMat.makeCompressed();
	std::cout << spMat << "\n";

	int n = spMat.rows() * spMat.cols(), nnz = 0;

	std::vector<int> pcols(n+1), irow;
	std::vector<std::complex<double>> values;
	pcols[0] = 0;

	std::vector<Eigen::Triplet<double>> triplets;
	for (int k = 0; k < spMat.outerSize(); ++k) {
		for (SparseComplexMatrix::InnerIterator it(spMat, k); it; ++it) {
			if (it.value() != std::complex<double>(0.0, 0.0)) {
				std::cout << "(" << it.row() << ", " << it.col() << ") = "
						<< it.value() << "\n";

				irow.push_back(it.row());
				values.push_back(it.value());
				triplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), std::real(it.value())));
				nnz++;
			}
		}
	}

	SparseRealMatrix spRealMat(spMat.rows(), spMat.cols());
	spRealMat.setFromTriplets(triplets.begin(), triplets.end());

	return std::real(0.0);
}

}
}




