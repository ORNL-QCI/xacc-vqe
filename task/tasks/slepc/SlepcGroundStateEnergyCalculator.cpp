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
	auto nTerms = inst.nInstructions();

	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	auto getBitStrForIdx = [&](std::uint64_t i) {
		std::stringstream s;
		for (int k = nQubits - 1; k >= 0; k--) s << ((i >> k) & 1);
		return s.str();
	};

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


	if (rank == 0) XACCInfo(
			"Building Matrix for SLEPc.");

	for (std::uint64_t myRow = Istart; myRow < Iend; myRow++) {

		auto rowBitStr = getBitStrForIdx(myRow);
		for (int i = 0; i < nTerms; i++) {

			std::shared_ptr<SpinInstruction> spinInst =
					std::dynamic_pointer_cast<SpinInstruction>(
							inst.getInstruction(i));
			auto newBitStrCoeff = spinInst->computeActionOnBra(rowBitStr);
			std::uint64_t k = std::stol(newBitStrCoeff.first, nullptr, 2);
			MatSetValue(A, myRow, k, newBitStrCoeff.second, ADD_VALUES);
		}
	}

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	if (rank == 0) XACCInfo(
			"Done building Matrix for SLEPc.");


	if (nQubits < 5) MatView(A, PETSC_VIEWER_STDOUT_WORLD);

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
	s << std::setprecision(12) << gsReal;
	if (rank == 0) XACCInfo("Lowest Eigenvalue = " + s.str());
	return std::real(gsReal);
}

}
}
