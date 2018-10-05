#include "XACC.hpp"
#include "VQEProgram.hpp"
#include "VQETask.hpp"
#include "VQEParameterGenerator.hpp"
#include "MPIProvider.hpp"
#include <boost/filesystem.hpp>

int main(int argc, char** argv) {

	// All our important stuff is in the xacc::vqe namespace
	using namespace xacc::vqe;

	// Initialize MPI, create a VQEProgram pointer
	std::shared_ptr<VQEProgram> program;

	// Add some command line options for this XACC app
	auto vqeOptions = std::make_shared<options_description>("XACC-VQE Options");
	vqeOptions->add_options()
				("vqe-program,f",value<std::string>(), "(required) The file containing "
							"XACC Kernels describing the Hamiltonian to simulate.")
				("vqe-task,t", value<std::string>(), "(required) The XACC-VQE Task to run.")
				("vqe-list-tasks,l","List available VQE Tasks.")
				("vqe-ansatz,a", value<std::string>(),"Provide the file name of the ansatz circuit.")
				("n-qubits,n",  value<std::string>(),"The number of qubits in the calculation")
				("n-electrons,e",  value<std::string>(),"The number of electrons in the calculation")
				("vqe-parameters,p",  value<std::string>(),"The initial parameters to seed VQE with, pass as string of comma separated parameters.")
				("vqe-energy-delta,d", value<std::string>(), "The change in energy to consider during classsical optimization.")
				("correct-readout-errors", "Correct qubit readout errors.")
				("qubit-map", "Provide a list of qubit indices as a comma-separated "
						"string to use in this computation. The 0th integer corresponds "
						"to the 0th logical qubit, etc.");

	xacc::addCommandLineOptions(vqeOptions);

	// We want to preempt the logger system with a global
	// predicate that only prints to rank == 0. So here we
	// tell the framework not to print during startup, but
	// to queue it so that when we set the predicate it will be dumped.
	std::vector<std::string> new_argv(argv, argv + argc);
//	new_argv.push_back("--queue-preamble");

	// Initialize the framework
	xacc::Initialize(new_argv);

	// Get the correct MPI Provider
	std::shared_ptr<MPIProvider> provider;
	std::shared_ptr<Communicator> world;
	if (xacc::hasService<MPIProvider>("boost-mpi")) {
		provider = xacc::getService<MPIProvider>("boost-mpi");
		provider->initialize(argc,argv);
		world = provider->getCommunicator();
		int rank = world->rank();
		xacc::setGlobalLoggerPredicate( [&]() { return rank == 0;});
		xacc::info("Using Boost MPI for distributed computations.");
	} else {
		provider = xacc::getService<MPIProvider>("no-mpi");
		provider->initialize(argc,argv);
		world = provider->getCommunicator();
		xacc::info("XACC-VQE Built without MPI Support.");
	}

	xacc::info("Number of Ranks = " + std::to_string(world->size()));
	if (!xacc::optionExists("accelerator")) {
		xacc::setAccelerator("vqe-dummy");
		// Set the default Accelerator to TNQVM
		if (xacc::hasAccelerator("tnqvm")) {
			xacc::setAccelerator("tnqvm");
		}
	}

	if (xacc::optionExists("vqe-list-tasks")) {
		xacc::info("");
		auto allTasks = xacc::getServices<VQETask>();
		for (auto a : allTasks) {
			xacc::info("VQE Task: " + a->name());
		}
		xacc::Finalize();
		exit(0);
	}

	// Users must specify a file containing VQE Hamiltonian kernels
	if (!xacc::optionExists("vqe-program")) {
		xacc::error("You must at least specify a kernel file to run this app.");
	}

	// Users must specify a task
	if (!xacc::optionExists("vqe-task")) {
		xacc::info("No VQE Task specified, defaulting to 'vqe-diagonalize'");
		xacc::setOption("vqe-task", "vqe-diagonalize");
	}

	// Get the user-specified Accelerator,
	// or TNQVM if none specified
	auto accelerator = xacc::getAccelerator();

	// Get the task to run
	auto task = xacc::getOption("vqe-task");

	// Read in the Hamiltonian kernel file
    auto fileName = xacc::getOption("vqe-program");
    if (!boost::filesystem::exists(fileName)) {
        xacc::error("Error: No such file - " + fileName);
    }
	std::ifstream moleculeKernelHpp(xacc::getOption("vqe-program"));
	std::string src((std::istreambuf_iterator<char>(moleculeKernelHpp)),
			std::istreambuf_iterator<char>());

	if (xacc::optionExists("vqe-ansatz")) {
        fileName = xacc::getOption("vqe-ansatz");
        if (!boost::filesystem::exists(fileName)) {
            xacc::error("Error: No such file - " + fileName);
        }
		std::ifstream spKernelHpp(xacc::getOption("vqe-ansatz"));
		std::string spsrc((std::istreambuf_iterator<char>(spKernelHpp)),
				std::istreambuf_iterator<char>());
		program = std::make_shared<VQEProgram>(accelerator, src, spsrc, world);
	} else {
		program = std::make_shared<VQEProgram>(accelerator, src, world);
	}

	program->build();

    auto buffer = accelerator->createBuffer("q",program->getNQubits());
    program->setGlobalBuffer(buffer);
    
	auto parameters = VQEParameterGenerator::generateParameters(program->getNParameters(), world);
	auto vqeTask = xacc::getService<VQETask>(task);
	vqeTask->setVQEProgram(program);

	VQETaskResult result = vqeTask->execute(parameters);

	//std::string msg = "Energy = ";
	//for (auto r : result.results) {
	//	std::stringstream ss;
	//	ss << std::setprecision(12) << r.second << " at (" << r.first.transpose() << ")";
	//	if (world->rank() == 0) xacc::info(msg + ss.str());
	//}

	xacc::Finalize();
}




