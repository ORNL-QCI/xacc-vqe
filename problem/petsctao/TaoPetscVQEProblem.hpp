#ifndef PROBLEM_TAOPETSCVQEPROBLEM_HPP_
#define PROBLEM_TAOPETSCVQEPROBLEM_HPP_

#include "VQEProblem.hpp"
#include <petsctao.h>

namespace xacc {

namespace vqe {

typedef struct {
  PetscInt  n;       /* dimension */
  PetscReal alpha;   /* condition parameter */
  PetscBool chained;
  VQEProblem* problem;
} AppCtx;

PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G,
		void *ptr) {
	PetscErrorCode ierr;
	AppCtx *user = (AppCtx *) ptr;
	double *x;

	/* Get pointers to vector data */
	ierr = VecGetArray(X, &x);

	std::cout << "HEY WORLD FROM TAO\n";
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

//		/* Allocate vectors for the solution and gradient */
		ierr = VecCreateSeq(PETSC_COMM_WORLD, user.n, &x);

		/* The TAO code begins here */

		/* Create TAO solver with desired solution method */
		ierr = TaoCreate(PETSC_COMM_WORLD, &tao);
		ierr = TaoSetType(tao, TAONM);

		/* Set solution vec and an initial guess */
		int *indices;
		VecSetValues(x, nParameters, indices, params.data());
		ierr = VecSet(x, zero);
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
