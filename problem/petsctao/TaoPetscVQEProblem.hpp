#ifndef PROBLEM_TAOPETSCVQEPROBLEM_HPP_
#define PROBLEM_TAOPETSCVQEPROBLEM_HPP_

#include "VQEProblem.hpp"
#include <petsctao.h>

namespace xacc {

namespace vqe {

typedef struct {
  PetscInt  n;
  PetscReal alpha;
  PetscBool chained;
  VQEProblem* problem;
  mpi::communicator comm;
} AppCtx;

PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G,
		void *ptr) {
	PetscErrorCode ierr;
	AppCtx *user = (AppCtx *) ptr;
	double *x;
	int localSize;

	/* Get pointers to vector data */
	ierr = VecGetArray(X, &x);
	ierr = VecGetLocalSize(X, &localSize);

	std::vector<double> allDataVec;
	std::vector<double> dataVec(x, x+localSize);

	boost::mpi::all_gather(user->comm, dataVec, &allDataVec);

//	if (user->comm.rank() == 0) std::cout << "Alldata vec:\n";
//	for (auto d : allDataVec) {
//		if (user->comm.rank() == 0) std::cout << d << " ";
//	}
//	if (user->comm.rank() == 0) std::cout << "\n";

	// Need to broadcast data pointer
	auto params = Eigen::Map<Eigen::VectorXd>(x, user->problem->getNParameters());
	*f = user->problem->value(params);

	/* Restore vectors */
	ierr = VecRestoreArray(X, &x);

	return 0;
}

class TaoPetscVQEProblem : public VQEProblem {

public:

	TaoPetscVQEProblem(std::istream& moleculeKernel,
			boost::mpi::communicator& c) :
			VQEProblem(moleculeKernel, c) {
	}

	TaoPetscVQEProblem(std::istream& moleculeKernel) :
			VQEProblem(moleculeKernel) {
	}

	virtual const Eigen::VectorXd minimize() {
		auto params = initializeParameters();

		PetscErrorCode ierr;
		PetscReal zero = 0.0;
		Vec x;
		Tao tao;
		PetscBool flg;
		PetscMPIInt size = comm.size(), rank = comm.rank();
		AppCtx user;

		/* Initialize TAO and PETSc */
		int argc = xacc::getArgc();
		char** argv = xacc::getArgv();
		ierr = PetscInitialize(&argc, &argv, (char*) 0, (char*) 0);

		/* Initialize problem parameters */
		user.n = nParameters;
		user.alpha = 99.0;
		user.chained = PETSC_FALSE;
		user.problem = this;
		user.comm = comm;

//		/* Allocate vectors for the solution and gradient */
//		ierr = VecCreateSeq(PETSC_COMM_WORLD, user.n, &x);
		ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, user.n, &x);

		/* The TAO code begins here */

		/* Create TAO solver with desired solution method */
		ierr = TaoCreate(PETSC_COMM_WORLD, &tao);
		ierr = TaoSetType(tao, TAONM);

		/* Set solution vec and an initial guess */
		int *indices = new int[nParameters];
		for (int i = 0; i < nParameters; i++) indices[i] = i;

		VecSetValues(x, nParameters, indices, params.data(), INSERT_VALUES);
		VecAssemblyBegin(x);
		VecAssemblyEnd(x);

		ierr = TaoSetInitialVector(tao, x);

		/* Set routines for function, gradient, hessian evaluation */
		ierr = TaoSetObjectiveAndGradientRoutine(tao, FormFunctionGradient,
				&user);

		/* Check for TAO command line options */
		ierr = TaoSetFromOptions(tao);

		/* SOLVE THE APPLICATION */
		ierr = TaoSolve(tao);

		ierr = TaoDestroy(&tao);
		ierr = VecDestroy(&x);

		return params;
	}
};
}
}

#endif
