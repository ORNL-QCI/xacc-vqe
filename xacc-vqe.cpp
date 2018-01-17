#include "XACC.hpp"
#include "ServiceRegistry.hpp"
#include "VQEProgram.hpp"
#include "VQETask.hpp"
#include "VQEParameterGenerator.hpp"

#include <boost/mpi.hpp>

int main(int argc, char** argv) {

	// All our important stuff is in the xacc::vqe namespace
	using namespace xacc::vqe;

	// Initialize MPI, create a VQEProgram pointer
	std::shared_ptr<VQEProgram> program;
	mpi::environment env(argc, argv);
	mpi::communicator world;

	xacc::setAccelerator("vqe-dummy");
	// Set the default Accelerator to TNQVM
	if (xacc::hasAccelerator("tnqvm")) {
		xacc::setAccelerator("tnqvm");
	}

	// Add some command line options for this XACC app

	auto vqeOptions = std::make_shared<options_description>("XACC-VQE Options");
	vqeOptions->add_options()
				("vqe-kernel-file,f",value<std::string>(), "(required) The file containing "
							"XACC Kernels describing a Hamiltonian.")
				("vqe-kernels-compiler,c", value<std::string>(), "")
				("vqe-task,t", value<std::string>(), "(required)")
				("vqe-list-tasks,l","List available VQE Tasks.")
				("vqe-state-prep-kernel,p", value<std::string>(),"(optional) Provide the file name of the state "
				"preparation circuit.")
				("vqe-state-prep-kernel-compiler,s",  value<std::string>(),"(optional) If not scaffold, provide the "
				"compiler to use in compiling the state prep circuit")
				("n-qubits,n",  value<std::string>(),"The number of qubits in the calculation")
				("n-electrons,e",  value<std::string>(),"The number of electrons in the calculation")
				("vqe-parameters",  value<std::string>(),"")
				("vqe-energy-delta", value<std::string>(), "");

	 xacc::addCommandLineOptions(vqeOptions);

	// Initialize the framework
	xacc::Initialize(argc, argv);

	if (xacc::optionExists("vqe-list-tasks")) {
		XACCInfo("");
		auto allTasks = xacc::ServiceRegistry::instance()->getServices<VQETask>();
		for (auto a : allTasks) {
			XACCInfo("VQE Task: " + a->name());
		}
		xacc::Finalize();
		exit(0);
	}

	// Users must specify a task
	if (!xacc::optionExists("vqe-task")) {
		XACCError("You must specify a vqe task to run.");
	}

	// Users must specify a file containing VQE Hamiltonian kernels
	if (!xacc::optionExists("vqe-kernel-file")) {
		XACCError("You must at least specify a kernel file to run this app.");
	}

	// Get the user-specified Accelerator,
	// or TNQVM if none specified
	auto accelerator = xacc::getAccelerator();

	// Get the task to run
	auto task = xacc::getOption("vqe-task");

	// Read in the Hamiltonian kernel file
	std::ifstream moleculeKernelHpp(xacc::getOption("vqe-kernel-file"));
	std::string src((std::istreambuf_iterator<char>(moleculeKernelHpp)),
			std::istreambuf_iterator<char>());

	if (xacc::optionExists("vqe-state-prep-kernel")) {
		std::ifstream spKernelHpp(xacc::getOption("vqe-state-prep-kernel"));
		std::string spsrc((std::istreambuf_iterator<char>(spKernelHpp)),
				std::istreambuf_iterator<char>());
		program = std::make_shared<VQEProgram>(accelerator, src, spsrc, world);
	} else {
		program = std::make_shared<VQEProgram>(accelerator, src, world);
	}

	program->build();

	auto parameters = VQEParameterGenerator::generateParameters(program->getNParameters(), world);
	auto vqeTask = xacc::ServiceRegistry::instance()->getService<VQETask>(task);
	vqeTask->setVQEProgram(program);

	VQETaskResult result = vqeTask->execute(parameters);

	std::string msg = (task == "compute-expectation-values" ? "Expectation Value = " : "Energy = ");
	for (auto r : result) {
		std::stringstream ss;
		ss << std::setprecision(12) << r.second << " at (" << r.first.transpose() << ")";
		if (world.rank() == 0) XACCInfo(msg + ss.str());
	}

	xacc::Finalize();
}




