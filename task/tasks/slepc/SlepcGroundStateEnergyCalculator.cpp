#include "SlepcGroundStateEnergyCalculator.hpp"
#include <petscmat.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <slepceps.h>
#include <unsupported/Eigen/KroneckerProduct>
#include "SpinInstruction.hpp"

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <memory>

using IndexPair = std::pair<std::uint64_t, std::uint64_t>;

namespace xacc {
namespace vqe {
double SlepcGroundStateEnergyCalculator::computeGroundStateEnergy(
		CompositeSpinInstruction& inst, const int nQubits) {
	boost::mpi::communicator world;
	int rank = world.rank();
	int nRanks = world.size();
	std::complex<double> gsReal;
	static char help[] = "";
	std::vector<std::string> argvVec;
	std::vector<char*> cstrs;
	argvVec.insert(argvVec.begin(), "appExec");
	for (auto& s : argvVec) {
		cstrs.push_back(&s.front());
	}

	int argc = argvVec.size();
	auto argv = cstrs.data();

	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	// Generate all bit strings
	std::vector<std::string> bitStrings(dim);
#pragma omp parallel for
	for (std::uint64_t j = 0; j < dim; j++) {

		std::stringstream s;
		for (int k = nQubits - 1; k >= 0; k--) {
			s << ((j >> k) & 1);
		}
		bitStrings[j] = s.str();
	}

	SlepcInitialize(&argc, &argv, (char*) 0, help);

	Mat A;
	EPS eps;
	EPSType type;
	PetscInt Istart, Iend;

	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, dim, dim);
	MatSetFromOptions(A);
	MatSetUp(A);

	MatGetOwnershipRange(A, &Istart, &Iend);

	auto nTerms = inst.nInstructions();

	// Get Identity coefficient
	std::complex<double> identityCoeff(0.0, 0.0);
	for (int i = 0; i < nTerms; i++) {
		std::shared_ptr<SpinInstruction> spinInst = std::dynamic_pointer_cast<
				SpinInstruction>(inst.getInstruction(i));
		if (spinInst->isIdentity()) {
			identityCoeff = spinInst->coefficient;
			break;
		}
	}

	if (rank == 0) XACCInfo(
			"Building Matrix for SLEPc.");
	for (std::uint64_t myRow = Istart; myRow < Iend; myRow++) {

		if (identityCoeff != std::complex<double>(0.0, 0.0)) {
			MatSetValue(A, myRow, myRow, identityCoeff, ADD_VALUES);
		}

		for (int i = 0; i < nTerms; i++) {

			std::shared_ptr<SpinInstruction> spinInst =
					std::dynamic_pointer_cast<SpinInstruction>(
							inst.getInstruction(i));
			std::pair<std::string, std::complex<double>> newBitStrCoeff;

			if (!spinInst->isIdentity()) {
				if (spinInst->isDiagonal()) {
					newBitStrCoeff = spinInst->computeActionOnBra(
							bitStrings[myRow]);
					MatSetValue(A, myRow, myRow, newBitStrCoeff.second,
							ADD_VALUES);
				} else {
					newBitStrCoeff = spinInst->computeActionOnBra(
							bitStrings[myRow]);
					std::uint64_t k = std::stol(newBitStrCoeff.first, nullptr,
							2);
					MatSetValue(A, myRow, k, newBitStrCoeff.second, ADD_VALUES);
				}
			}
		}
	}

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	if (rank == 0) XACCInfo(
			"Done building Matrix for SLEPc.");


//	MatView(A, PETSC_VIEWER_STDOUT_WORLD);

	EPSCreate(PETSC_COMM_WORLD, &eps);
	EPSSetOperators(eps, A, NULL);
	EPSSetProblemType(eps, EPS_HEP);
	EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
	EPSSetType(eps, EPSLANCZOS);
	if (rank == 0)
		XACCInfo("Starting SLEPc EPS Solve");
	EPSSolve(eps);
	if (rank == 0)
		XACCInfo("Done with SLEPc EPS Solve");
	EPSGetEigenpair(eps, 0, &gsReal, NULL, NULL, NULL);
	EPSDestroy(&eps);
	MatDestroy(&A);
	SlepcFinalize();

	std::stringstream s;
	s << gsReal;
	XACCInfo("Lowest Eigenvalue = " + s.str());
	return std::real(gsReal);
}

}
}
