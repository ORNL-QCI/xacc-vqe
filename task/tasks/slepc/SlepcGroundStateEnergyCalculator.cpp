#include "SlepcGroundStateEnergyCalculator.hpp"
#include <petscmat.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <slepceps.h>

namespace xacc {
namespace vqe {
double SlepcGroundStateEnergyCalculator::computeGroundStateEnergy(
		CompositeSpinInstruction& inst, const int nQubits) {
	boost::mpi::communicator world;
	int rank = world.rank();
	int nRanks = world.size();
	double gsReal;
	static char help[] = "";
	std::vector<std::string> argvVec;
	std::vector<char*> cstrs;
	argvVec.insert(argvVec.begin(), "appExec");
	for (auto& s : argvVec) {
		cstrs.push_back(&s.front());
	}

	int argc = argvVec.size();
	auto argv = cstrs.data();

	using SparseComplexMatrix = Eigen::SparseMatrix<std::complex<double>>;
	SparseComplexMatrix spMat = inst.toSparseMatrix(nQubits).pruned();
	spMat.makeCompressed();

	Mat A;
	EPS eps;
	EPSType type;
	PetscInt n = spMat.rows(), Istart, Iend;

	SlepcInitialize(&argc, &argv, (char*) 0, help);

	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
	MatSetFromOptions(A);
	MatSetUp(A);

	MatGetOwnershipRange(A, &Istart, &Iend);
	for (int k = 0; k < spMat.outerSize(); ++k) {
		for (SparseComplexMatrix::InnerIterator it(spMat, k); it; ++it) {
			int r = it.row();
			int col = it.col();
			double val = std::real(it.value());
			// Only add our rank's values
			if (r >= Istart && r < Iend) {
				MatSetValue(A, r, col, val, INSERT_VALUES);
			}
		}
	}

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

//	MatView(A, PETSC_VIEWER_STDOUT_WORLD);
	EPSCreate(PETSC_COMM_WORLD, &eps);
	EPSSetOperators(eps, A, NULL);
	EPSSetProblemType(eps, EPS_HEP);
	EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
	EPSSetType(eps, EPSLANCZOS);
	EPSSolve(eps);
	EPSGetEigenpair(eps, 0, &gsReal, NULL, NULL, NULL);
	EPSDestroy(&eps);
	MatDestroy(&A);
	SlepcFinalize();

	return std::real(gsReal);
}

}
}

