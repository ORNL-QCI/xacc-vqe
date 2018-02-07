#include <petscmat.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <petsctao.h>

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <memory>
#include "TaoVQEBackend.hpp"
#include "MPIProvider.hpp"

using IndexPair = std::pair<std::uint64_t, std::uint64_t>;

namespace xacc {
namespace vqe {

typedef struct {
  int nParameters;
  std::shared_ptr<ComputeEnergyVQETask> computeTask;
  double currentEnergy = 0.0;
} AppCtx;

PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G,
		void *ptr) {
	AppCtx *user = (AppCtx *) ptr;
	const double* x;

	/* Get pointers to vector data */
	VecGetArrayRead(X, &x);

	// Need to broadcast data pointer
	auto params = Eigen::Map<const Eigen::VectorXd>(x, user->nParameters);
	auto e = user->computeTask->execute(params).energy;
	*f = e;
	user->currentEnergy = e;

	/* Restore vectors */
	VecRestoreArrayRead(X, &x);

	return 0;
}

const VQETaskResult TaoVQEBackend::minimize(Eigen::VectorXd parameters) {

	static char help[] = "";
	std::vector<std::string> argvVec;
	std::vector<char*> cstrs;
	argvVec.insert(argvVec.begin(), "appExec");
	for (auto& s : argvVec) {
		cstrs.push_back(&s.front());
	}

	int argc = argvVec.size();
	auto argv = cstrs.data();

	auto nParameters = program->getNParameters();

	std::cout << "HI: " << parameters << "\n";
	PetscErrorCode ierr;
	Vec x;
	Tao tao;
	PetscBool flg;
	AppCtx user;

	computeTask = std::make_shared<ComputeEnergyVQETask>(program);

	/* Initialize TAO and PETSc */
	ierr = PetscInitialize(&argc, &argv, (char*) 0, (char*) 0);

	/* Initialize problem parameters */
	user.nParameters = nParameters;
	user.computeTask = computeTask;

	ierr = VecCreateSeq(PETSC_COMM_WORLD, nParameters, &x);

	ierr = TaoCreate(PETSC_COMM_WORLD, &tao);
	ierr = TaoSetType(tao, TAONM);

	/* Set solution vec and an initial guess */
	PetscInt *indices = new PetscInt[nParameters];
	for (int i = 0; i < nParameters; i++) indices[i] = i;

//	std::vector<std::complex<double>> complexParams;
//	for (int i = 0; i < parameters.cols(); i++) {
//		complexParams.push_back(std::complex<double>(parameters(i),0));
//	}

	VecSetValues(x, nParameters, indices, parameters.data(), INSERT_VALUES);
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);

	ierr = TaoSetInitialVector(tao, x);

	/* Set routines for function, gradient, hessian evaluation */
	ierr = TaoSetObjectiveAndGradientRoutine(tao, FormFunctionGradient,
			&user);

	ierr = TaoSolve(tao);

	ierr = TaoDestroy(&tao);
	ierr = VecDestroy(&x);

	PetscFinalize();
	VQETaskResult result;
	result.energy = user.currentEnergy;

	return result;
}

}
}
