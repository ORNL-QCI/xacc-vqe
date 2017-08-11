#include "XACC.hpp"
#include "VQEProblem.hpp"
#include "solver/neldermeadsolver.h"

int main(int argc, char** argv) {

	using namespace xacc::vqe;

	xacc::Initialize(argc, argv);

	xacc::addCommandLineOption("kernel-directory", "The directory "
			"containing *.hpp files, each containing an XACC Kernel describing "
			"a molecular fermionic Hamiltonian.");

	std::ifstream moleculeKernelHpp(
			"../vqe/scripts/H2_sto-3g_singlet_H2_Molecule.hpp");

	VQEProblem<double> problem(moleculeKernelHpp);

	auto params = problem.initializeParameters();

	std::cout << "PARAMS VEC:\n" << params << "\n";

	cppoptlib::NelderMeadSolver<VQEProblem<double>> solver;

	solver.minimize(problem, params);

	std::stringstream ss;
	ss << params.transpose();
	XACCInfo(std::string(42, '-'));
	XACCInfo("Final VQE_Params: (" + ss.str() + ")");
	problem(params);

	xacc::Finalize();

}

