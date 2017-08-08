#include "XACC.hpp"
#include "VQEProblem.hpp"
#include "solver/neldermeadsolver.h"

int main(int argc, char** argv) {

	using namespace xacc::vqe;

	xacc::Initialize(argc, argv);

	std::ifstream moleculeKernelHpp(
			"../vqe/scripts/H2_sto-3g_singlet_H2_Molecule.hpp");

	VQEProblem<double> problem(moleculeKernelHpp);

	auto params = problem.initializeParameters();

	std::cout << "PARAMS VEC:\n" << params << "\n";

	cppoptlib::NelderMeadSolver<VQEProblem<double>> solver;

	solver.minimize(problem, params);

	std::cout << std::string(42, '-') << std::endl;
	std::cout << "   argmin: " << params.transpose() << std::endl;
	std::cout << "   f in argmin: " << problem(params) << std::endl;

	xacc::Finalize();

}

