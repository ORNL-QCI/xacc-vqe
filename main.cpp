#include "XACC.hpp"
#include "VQEProblem.hpp"
#include "solver/neldermeadsolver.h"

int main(int argc, char** argv) {

	using namespace xacc::vqe;

	xacc::addCommandLineOption("kernel-directory", "The directory "
			"containing *.hpp files, each containing an XACC Kernel describing "
			"a molecular fermionic Hamiltonian.");
	xacc::addCommandLineOption("kernel-file", "The file "
			"containing an XACC Kernel describing "
			"a molecular fermionic Hamiltonian.");

	xacc::Initialize(argc, argv);



	auto options = xacc::RuntimeOptions::instance();

	if (options->exists("kernel-directory")) {
		std::vector<std::ifstream> fileStreams;
		std::vector<double> moleculeParams;
		std::string outFileName;
		boost::filesystem::path filePath((*options)["kernel-directory"]);
		boost::filesystem::directory_iterator end_itr;
		for (boost::filesystem::directory_iterator itr(filePath);
				itr != end_itr; ++itr) {
			auto p = itr->path();
			if (p.extension() == ".hpp") {
				outFileName = p.parent_path().filename().string();
				fileStreams.push_back(std::ifstream(p.string()));
				std::vector<std::string> splitUnderscore;
				boost::split(splitUnderscore, p.string(), boost::is_any_of("_"));
				// should be the last one in the vector
				auto param = splitUnderscore[splitUnderscore.size()-1];
				boost::replace_all(param, ".hpp", "");
				moleculeParams.push_back(std::stod(param));
			}
		}

		// We have all files, now lets run VQE on all of them
		std::ofstream outFile(outFileName+".csv");
		for (int i = 0; i < fileStreams.size(); i++) {
			VQEProblem<double> problem(fileStreams[i]);
			auto params = problem.initializeParameters();
			cppoptlib::NelderMeadSolver<VQEProblem<double>> solver;
			solver.minimize(problem, params);
			auto finalValue = problem(params);
			outFile << moleculeParams[i] << ", " << finalValue << "\n";
		}

		outFile.close();

	} else {

		if (!options->exists("kernel-file")) {
			XACCError("You must at least specify a kernel file to run this app.");
		}

		std::ifstream moleculeKernelHpp((*options)["kernel-file"]);

		VQEProblem<double> problem(moleculeKernelHpp);

		auto params = problem.initializeParameters();

		cppoptlib::NelderMeadSolver<VQEProblem<double>> solver;

		solver.minimize(problem, params);

		std::stringstream ss;
		ss << params.transpose();
		XACCInfo(std::string(42, '-'));
		XACCInfo("Final VQE_Params: (" + ss.str() + ")");
		problem(params);
	}

	xacc::Finalize();

}

