#include "XACC.hpp"
#include "VQEProblem.hpp"
#include "solver/neldermeadsolver.h"

using ExecutionPair = std::pair<double, std::string>;
struct ExecutionPairComparison {
	constexpr bool operator()(const ExecutionPair& lhs,
			const ExecutionPair& rhs) const {
		return rhs.first < lhs.first;
	}
};

int main(int argc, char** argv) {

	// All our important stuff is in the xacc::vqe namespace
	using namespace xacc::vqe;

	// Add some command line options for XACC VQE
	xacc::addCommandLineOptions("XACC VQE", std::map<std::string, std::string>{
		{"kernel-directory", "The directory containing *.hpp files, each "
				"containing an XACC Kernel describing a molecular "
				"fermionic Hamiltonian."},
		{"kernel-file", "The file containing an XACC Kernel describing "
			"a molecular fermionic Hamiltonian."},
		{"vqe-energy-delta", "Specify the change in the energy condition for convergence of Nelder-Mead minimization."},
		{"vqe-iterations", "The number of iterations before stoping Nelder-Mead minimization."}
	});

	// Initialize the Framework
	xacc::Initialize(argc, argv);

	// We can either do a parameterized execution, ie,
	// run a bunch of calculations, one for each file in
	// the provided directory,
	// OR
	// Run a calculation on oone file.
	if (xacc::optionExists("kernel-directory")) {
		std::priority_queue<ExecutionPair, std::vector<ExecutionPair>,
				ExecutionPairComparison> executions;
		std::string outFileName;
		boost::filesystem::path filePath(xacc::getOption("kernel-directory"));
		boost::filesystem::directory_iterator end_itr;
		for (boost::filesystem::directory_iterator itr(filePath);
				itr != end_itr; ++itr) {
			auto p = itr->path();
			if (p.extension() == ".hpp") {
				outFileName = p.parent_path().filename().string();
				std::vector<std::string> splitUnderscore;
				boost::split(splitUnderscore, p.string(), boost::is_any_of("_"));
				// should be the last one in the vector
				auto param = splitUnderscore[splitUnderscore.size()-1];
				boost::replace_all(param, ".hpp", "");
				executions.push({std::stod(param), p.string()});
			}
		}

		// We have all files, now lets run VQE on all of them
		std::ofstream outFile(outFileName+".csv");
		Eigen::VectorXd params(2);
		params << 1.69, -.433;
		bool executedOnce = false;
		while(!executions.empty()) {
			std::ifstream stream(executions.top().second);
			VQEProblem<double> problem(stream);
//			if (!executedOnce) {
//				// We do this so that at first we generate
//				// random parameters, but then we just keep
//				// reusing the last iterations parameters...
//				params = problem.initializeParameters();
//				executedOnce = true;
//			}
			cppoptlib::NelderMeadSolver<VQEProblem<double>> solver;
			solver.setStopCriteria(VQEProblem<double>::getConvergenceCriteria());
			std::cout << "Starting Parameters: " << params.transpose() << "\n";
			solver.minimize(problem, params);
			auto finalValue = problem.currentEnergy;
			std::cout << "Converged to " << finalValue << " with " << params.transpose() << "\n";
			outFile << executions.top().first << ", " << finalValue << "\n";
			outFile.flush();
			executions.pop();
			stream.close();
		}

		outFile.close();

	} else {

		if (!xacc::optionExists("kernel-file")) {
			XACCError("You must at least specify a kernel file to run this app.");
		}

		std::ifstream moleculeKernelHpp(xacc::getOption("kernel-file"));

		VQEProblem<double> problem(moleculeKernelHpp);

		auto params = problem.initializeParameters();

		params << 1.69, -.433;

		cppoptlib::NelderMeadSolver<VQEProblem<double>> solver;
		solver.setStopCriteria(VQEProblem<double>::getConvergenceCriteria());

		solver.minimize(problem, params);

		std::stringstream ss;
		ss << params.transpose();
		XACCInfo(std::string(42, '-'));
		XACCInfo("Final VQE_Params: (" + ss.str() + ")");
		problem(params);
	}

	xacc::Finalize();

}

