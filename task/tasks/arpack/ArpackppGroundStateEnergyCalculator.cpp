#include "ArpackppGroundStateEnergyCalculator.hpp"
#include "arlsmat.h"
#include "arlssym.h"
#include "arssym.h"

namespace xacc {
namespace vqe {
double ArpackppGroundStateEnergyCalculator::computeGroundStateEnergy(
		CompositeSpinInstruction& inst, const int nQubits) {
	boost::mpi::communicator world;

//	using SparseMatrix = Eigen::SparseMatrix<std::complex<double>>;
//	SparseMatrix spMat = inst.toSparseMatrix(nQubits);
//	spMat.makeCompressed();
//	std::cout << spMat << "\n";
//
//	int n = spMat.rows() * spMat.cols(), nnz = 0;
//
//	std::vector<int> pcols(n+1), irow;
//	std::vector<std::complex<double>> values;
//	pcols[0] = 0;
//
//	for (int k = 0; k < spMat.outerSize(); ++k) {
//		for (SparseMatrix::InnerIterator it(spMat, k); it; ++it) {
//			if (it.value() != std::complex<double>(0.0, 0.0)) {
//				std::cout << "(" << it.row() << ", " << it.col() << ") = "
//						<< it.value() << "\n";
//
//				irow.push_back(it.row());
//				values.push_back(it.value());
//
//				nnz++;
//			}
//		}
//	}
//
//	pcols[n] = nnz;
//
//	std::cout << "NONzeros: " << nnz << ", " << spMat.rows()*spMat.cols() << "\n";
//
//	char uplo = 'U';
//	ARluSymMatrix<std::complex<double>> A(spMat.cols(), nnz, values.data(), irow.data(),
//			pcols.data(), uplo);
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
	spRealMat.makeCompressed();

	MyMatrix<double> A(spRealMat);

	ARSymStdEig<double, MyMatrix<double>> eig(spRealMat.cols(), 4L, &A,
			&MyMatrix<double>::MultMv, "SM");
	eig.FindEigenvectors();

	n = eig.GetN();
	int nconv = eig.ConvergedEigenvalues();
	int mode = eig.GetMode();

	std::cout << std::endl << std::endl
			<< "Testing ARPACK++ class ARSymStdEig \n";
	std::cout << "Real symmetric eigenvalue problem: A*x - lambda*x"
			<< std::endl;
	switch (mode) {
	case 1:
		std::cout << "Regular mode" << std::endl << std::endl;
		break;
	case 3:
		std::cout << "Shift and invert mode" << std::endl << std::endl;
	}

	std::cout << "Dimension of the system            : " << n << std::endl;
	std::cout << "Number of 'requested' eigenvalues  : " << eig.GetNev()
			<< std::endl;
	std::cout << "Number of 'converged' eigenvalues  : " << nconv << std::endl;
	std::cout << "Number of Arnoldi vectors generated: " << eig.GetNcv()
			<< std::endl;
	std::cout << "Number of iterations taken         : " << eig.GetIter()
			<< std::endl;
	std::cout << std::endl;

	   for (int i=0; i<nconv; i++) {
			std::cout << "HELLO: " << eig.Eigenvalue(i) << "\n";
	    }

		Eigen::SelfAdjointEigenSolver<SparseRealMatrix> es(spRealMat);
		auto eigenvalues = es.eigenvalues();
		std::stringstream ss;
		ss << eigenvalues.transpose();
		if (world.rank() == 0)
			XACCInfo("HamiltonianEigenvalues:\n" + ss.str());
//	ARluCompStdEig<double> prob(1, A, "SM");
	return std::real(0.0);
}

}
}




