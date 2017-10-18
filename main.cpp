#include "XACC.hpp"
#include "VQEProblem.hpp"
#include "solver/neldermeadsolver.h"
#include "solver/conjugatedgradientdescentsolver.h"
#include "solver/gradientdescentsolver.h"

#include <boost/mpi.hpp>

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

	mpi::environment env(argc, argv);
	mpi::communicator world;

	// Set the default Accelerator to TNQVM, and
	// default number of electrons to 2
	xacc::setAccelerator("tnqvm");
	xacc::setOption("n-electrons", "2");

	// Add some command line options for XACC VQE
	xacc::addCommandLineOptions("XACC VQE", std::map<std::string, std::string>{
		{"vqe-kernel-directory", "The directory containing *.hpp files, each "
				"containing an XACC Kernel describing a molecular "
				"fermionic Hamiltonian."},
		{"vqe-kernel-file", "The file containing an XACC Kernel describing "
			"a molecular fermionic Hamiltonian."},
		{"vqe-energy-delta", "Specify the change in the energy condition for "
				"convergence of Nelder-Mead minimization."},
		{"vqe-iterations", "The number of iterations before stoping "
				"Nelder-Mead minimization."},
		{"vqe-print-scaffold-source", "Print the source code in the Scaffold "
				"language to the provided file name."},
		{"vqe-exit-after-scaffold", "Print Scaffold source code to the "
				"provided file name and then quit. Set this arg to be y"},
		{"vqe-state-prep-kernel", "Provide the file name of the state "
				"preparation circuit."},
		{"vqe-state-prep-kernel-compiler", "If not scaffold, provide the "
				"compiler to use in compiling the state prep circuit."},
		{"n-qubits", "The number of qubits used in this calculation."},
		{"vqe-print-stats", "Print the number of qubits, variational parameters, and Hamiltonian terms."}
	});

	std::cout << "PID: " << world.rank() << ", " << getpid() << "\n";
	// Initialize the Framework
	xacc::Initialize(argc, argv);

	// We can either do a parameterized execution, ie,
	// run a bunch of calculations, one for each file in
	// the provided directory,
	// OR
	// Run a calculation on oone file.
	if (xacc::optionExists("vqe-kernel-directory")) {
		std::priority_queue<ExecutionPair, std::vector<ExecutionPair>,
				ExecutionPairComparison> executions;
		std::string outFileName;
		boost::filesystem::path filePath(xacc::getOption("vqe-kernel-directory"));
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
		bool executedOnce = false;
		while(!executions.empty()) {
			std::ifstream stream(executions.top().second);
			VQEProblem problem(stream, world);
			params = problem.initializeParameters();
			cppoptlib::NelderMeadSolver<VQEProblem> solver;
			solver.setStopCriteria(VQEProblem::getConvergenceCriteria());
			solver.minimize(problem, params);
			auto finalValue = problem.currentEnergy;
			outFile << executions.top().first << ", " << finalValue << "\n";
			outFile.flush();
			executions.pop();
			stream.close();
			if (world.rank() == 0) XACCInfo("Number of Hamiltonian Terms = " + std::to_string(problem.kernels.size()));
			if (world.rank() == 0) XACCInfo("Total QPU Calls = " + std::to_string(problem.totalQpuCalls));
		}

		outFile.close();

	} else {

		if (!xacc::optionExists("vqe-kernel-file")) {
			XACCError("You must at least specify a kernel file to run this app.");
		}

		std::ifstream moleculeKernelHpp(xacc::getOption("vqe-kernel-file"));
		VQEProblem problem(moleculeKernelHpp, world);

		Eigen::VectorXd params = problem.initializeParameters();

		cppoptlib::NelderMeadSolver<VQEProblem> solver;
		solver.setStopCriteria(VQEProblem::getConvergenceCriteria());
		solver.minimize(problem, params);

		std::stringstream ss;
		ss << params.transpose();
		if (world.rank() == 0)XACCInfo(std::string(42, '-'));
		if (world.rank() == 0)XACCInfo("Final VQE_Params: (" + ss.str() + ")");
		problem(params);

		if (world.rank() == 0)XACCInfo("Number of Hamiltonian Terms = " + std::to_string(problem.kernels.size()));
		if (world.rank() == 0)XACCInfo("Total QPU Calls = " + std::to_string(problem.totalQpuCalls));
		if (world.rank() == 0)XACCInfo("Total VQE Steps = " + std::to_string(problem.totalQpuCalls/problem.kernels.size()));


	}

	xacc::Finalize();
}

