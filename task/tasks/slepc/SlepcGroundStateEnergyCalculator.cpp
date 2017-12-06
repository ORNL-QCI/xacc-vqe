#include "SlepcGroundStateEnergyCalculator.hpp"
#include <petscmat.h>
#include <petscsys.h>
#include <petscviewer.h>

namespace xacc {
namespace vqe {
double SlepcGroundStateEnergyCalculator::computeGroundStateEnergy(
		CompositeSpinInstruction& inst, const int nQubits) {
	boost::mpi::communicator world;
	int rank = world.rank();
	int nRanks = world.size();

	using SparseComplexMatrix = Eigen::SparseMatrix<std::complex<double>>;
	using SparseRealMatrix = Eigen::SparseMatrix<double>;
	SparseComplexMatrix spMat = inst.toSparseMatrix(nQubits).pruned();
	spMat.makeCompressed();
	if (rank == 0) std::cout << spMat << "\n";

	int myRowStart = (rank) * spMat.rows() / nRanks;
	int myRowEnd = (rank + 1) * spMat.rows() / nRanks;
	int myColStart = myRowStart;
	int myColEnd = myRowEnd;

	int localNRows = myRowEnd - myRowStart;
	int localNCols = myColEnd - myColStart;

	std::cout << "HELLO WORLD: " << myRowStart << ", " << myRowEnd << ", " << myColStart << ", " << myColEnd << "\n";
	std::vector<int> myRowIndices(localNRows), myColIndices(spMat.cols());
	std::iota(myRowIndices.begin(), myRowIndices.end(), myRowStart);
	std::iota(myColIndices.begin(), myColIndices.end(), 0);

	std::map<int, int> globalToLocalRowIdx, localToGlobalRowIdx, globalToLocalColIdx;
	int j = 0;
	for (int i = myRowStart; i < myRowEnd; i++) {
		globalToLocalRowIdx[i] = j;
		localToGlobalRowIdx[j] = i;
		j++;
	}

	j = 0;
	for (int i = myColStart; i < myColEnd; i++) {
		globalToLocalColIdx[i] = j;
		j++;
	}

	std::vector<int> idxm, idxn;
	std::vector<double> values(localNRows*localNCols);
	for (int k = 0; k < spMat.outerSize(); ++k) {
		for (SparseComplexMatrix::InnerIterator it(spMat, k); it; ++it) {

			// Only add our rank's values
			if (it.row() >= myRowStart && it.row() < myRowEnd) {
				if (rank == 0) std::cout << "(" << it.row() << ", " << it.col() << ") = "
						<< it.value() << "\n";

				idxm.push_back(it.row());
				idxn.push_back(it.col());
				int i = globalToLocalRowIdx[it.row()];
				int j = globalToLocalColIdx[it.col()];

				std::cout << "Value at " << i << ", " << j << " = " << std::real(it.value()) << "\n";
				values[i*localNCols + j] = std::real(it.value());
			}
		}
	}

	for(int i = 0; i < nRanks; i++) {
	    world.barrier();
	    if (i == rank) {
	        std::cout << rank << " = " << localNRows << ", " << localNCols << "\n";
	    }
	}

	SparseComplexMatrix diag = spMat.block(
			rank * ((rank == nRanks - 1) ? localNRows - 1 : localNRows),
			rank * ((rank == nRanks - 1) ? localNCols - 1 : localNCols),
			localNRows, localNCols);

	for(int i = 0; i < nRanks; i++) {
	    world.barrier();
	    if (i == rank) {
	        std::cout << rank << " = \n" << diag << "\n";
	    }
	}

	std::vector<int> counter(localNRows), totalCounter(localNRows), offDiagNonZeros(localNRows);
	auto inner = diag.innerIndexPtr();
	for (int i = 0; i < diag.nonZeros(); i++) {
		counter[inner[i]]++;
	}

//	int* d_nnz = counter.data();

	auto totalInner = spMat.innerIndexPtr();
	for (int i = 0; i < spMat.nonZeros(); i++) {
		totalCounter[totalInner[i]]++;
	}

	for (int i = 0; i < localNRows; i++) {
		offDiagNonZeros[i] = totalCounter[i] - counter[i];
	}

//	int* od_nnz = offDiagNonZeros.data();

	auto argc = xacc::getArgc();
	auto argv = xacc::getArgv();

	static char help[] = "";

	PetscInitialize(&argc, &argv, 0, help);

	Mat A;
//	MatCreateAIJ(PETSC_COMM_WORLD, localNRows, localNCols, spMat.rows(),
//			spMat.cols(), 0, d_nnz, 2, 0, &A);


	MatCreateAIJ(PETSC_COMM_WORLD, localNRows, localNCols, spMat.rows(),
			spMat.cols(), 0, NULL, 0, NULL, &A);

//	MatSetValues(A, localNRows, idxm.data(), localNCols, idxn.data(),
//			values.data(), INSERT_VALUES);

	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

//	MatView(A, PETSC_VIEWER_STDOUT_WORLD);

	MatDestroy(&A);
	PetscFinalize();

	return std::real(0.0);
}

}
}




