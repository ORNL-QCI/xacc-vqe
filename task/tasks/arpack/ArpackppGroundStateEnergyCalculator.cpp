#include "ArpackppGroundStateEnergyCalculator.hpp"

#include "arlnsmat.h"
#include "arlscomp.h"

namespace xacc {
namespace vqe {
double ArpackppGroundStateEnergyCalculator::computeGroundStateEnergy(
		CompositeSpinInstruction& inst, const int nQubits) {
	boost::mpi::communicator world;

	using SparseMatrix = Eigen::SparseMatrix<std::complex<double>>;
	SparseMatrix spMat = inst.toSparseMatrix(nQubits);
	spMat.makeCompressed();
	std::cout << spMat << "\n";

	std::vector<int> irow, pcol;
	int nnz = 0;
	std::vector<std::complex<double>> values;
	for (int k = 0; k < spMat.outerSize(); ++k) {
		for (SparseMatrix::InnerIterator it(spMat, k); it; ++it) {
			if (it.value() != std::complex<double>(0.0, 0.0)) {
				std::cout << "(" << it.row() << ", " << it.col() << ") = "
						<< it.value() << "\n";
				irow.push_back(it.row());
				pcol.push_back(it.col());
				values.push_back(it.value());
				nnz++;
			}
		}
	}

	std::cout << "NONzeros: " << nnz << ", " << spMat.rows()*spMat.cols() << "\n";

	ARluNonSymMatrix<std::complex<double>, double> A(
			spMat.rows() * spMat.cols(), nnz, values.data(), irow.data(),
			pcol.data());

	return std::real(0.0);
}

}
}




