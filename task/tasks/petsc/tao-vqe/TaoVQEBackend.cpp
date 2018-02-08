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
  Eigen::VectorXd observations;
  std::vector<Eigen::VectorXd> angles;
} AppCtx;

PetscErrorCode nelderMeadFunction(Tao tao, Vec X, PetscReal *f, Vec G,
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

PetscErrorCode poundersFunction(Tao tao, Vec X, Vec F, void * ptr) {
	AppCtx *user = (AppCtx *) ptr;
	PetscInt i;
	PetscReal* f, *x;

	VecGetArray(X, &x);
	VecGetArray(F, &f);

	auto params = Eigen::Map<const Eigen::VectorXd>(x, user->nParameters);

	auto y = user->observations;
	auto thetas = user->angles;
	auto e = user->computeTask->execute(params).energy;
	for (i = 0; i < y.rows(); i++) {
		f[i] = 0.5 * (y(i) - e);
		std::cout << "Residual " << i << ": " << f[i] << "\n";
	}
	VecRestoreArray(X, &x);
	VecRestoreArray(F, &f);
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
	PetscInitialize(&argc, &argv, (char*) 0, (char*) 0);

	auto nParameters = program->getNParameters();

	Vec x;
	Tao tao;
	PetscBool flg;
	AppCtx user;

	computeTask = std::make_shared<ComputeEnergyVQETask>(program);
	user.nParameters = nParameters;
	user.computeTask = computeTask;
	VecCreateSeq(PETSC_COMM_WORLD, nParameters, &x);
	TaoCreate(PETSC_COMM_WORLD, &tao);
	std::string t = "nm";
	if (xacc::optionExists("tao-type")) {
		t = xacc::getOption("tao-type");
	}
	TaoSetType(tao, t == "nm" ? TAONM : TAOPOUNDERS);

	/* Set solution vec and an initial guess */
	PetscInt *indices = new PetscInt[nParameters];
	for (int i = 0; i < nParameters; i++)
		indices[i] = i;

	// Set the Parameter vector
	VecSetValues(x, nParameters, indices, parameters.data(), INSERT_VALUES);
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);

	TaoSetInitialVector(tao, x);

	if (t == "pounders") {
		Vec F;

		if (!xacc::optionExists("observations")) {
			xacc::error(
					"You must provide a file containing previous observations for Pounders.");
		}
		auto obsFileName = xacc::getOption("observations");
		std::ifstream obsFile(obsFileName);
		std::string contents((std::istreambuf_iterator<char>(obsFile)),
				std::istreambuf_iterator<char>());
		std::vector<std::string> split;
		boost::split(split, contents, boost::is_any_of("\n"));
		Eigen::VectorXd energies(split.size()-1);
		std::vector<Eigen::VectorXd> angles;

		VecCreateSeq(PETSC_COMM_WORLD, energies.rows(), &F);

		int counter = 0;
		for (auto line : split) {
			if (!line.empty()) {
				std::vector<std::string> split2;
				boost::split(split2, line, boost::is_any_of(","));
				int nParams = split2.size() - 1;
				Eigen::VectorXd ps(nParams);
				for (int i = 0; i < nParams; i++) {
					ps(i) = std::stod(split2[i]);
				}

				angles.push_back(ps);
				energies(counter) = std::stod(split2[nParams]);
				counter++;
			}
		}

		user.angles = angles;
		user.observations = energies;

		TaoSetSeparableObjectiveRoutine(tao, F, poundersFunction,
				(void*) &user);

	} else {
		/* Set routines for function, gradient, hessian evaluation */
		TaoSetObjectiveAndGradientRoutine(tao, nelderMeadFunction, &user);
	}

	TaoSolve(tao);

	TaoDestroy(&tao);
	VecDestroy(&x);
	PetscFinalize();

	VQETaskResult result;
	result.energy = user.currentEnergy;
	return result;
}

}
}
