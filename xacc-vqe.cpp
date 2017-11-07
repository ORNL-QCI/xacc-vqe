#include "XACC.hpp"
#include "VQEProgram.hpp"
#include "VQETask.hpp"
#include "VQEParameterGenerator.hpp"

#include <boost/mpi.hpp>

int main(int argc, char** argv) {

	// All our important stuff is in the xacc::vqe namespace
	using namespace xacc::vqe;

	//VQETaskFactory factory;
	std::shared_ptr<VQEProgram> program;
	mpi::environment env(argc, argv);
	mpi::communicator world;

	// Set the default Accelerator to TNQVM
	xacc::setAccelerator("tnqvm");

	xacc::addCommandLineOptions("XACC VQE", std::map<std::string, std::string>{
		{"vqe-kernel-file", "(required) The file containing an XACC Kernel describing "
			"a molecular fermionic Hamiltonian."},
		{"vqe-task", "(required)"},
		{"vqe-state-prep-kernel", "(optional) Provide the file name of the state "
				"preparation circuit."},
		{"vqe-state-prep-kernel-compiler", "(optional) If not scaffold, provide the "
				"compiler to use in compiling the state prep circuit."},
		{"n-qubits", "The number of qubits used in this calculation."},
		{"vqe-initial-parameters", "(optional)"}
	});

	xacc::Initialize(argc, argv);

	if (!xacc::optionExists("vqe-task")) {
		XACCError("You must specify a vqe task to run.");
	}

	if (!xacc::optionExists("vqe-kernel-file")) {
		XACCError("You must at least specify a kernel file to run this app.");
	}

	auto accelerator = xacc::getAccelerator();
	auto task = xacc::getOption("vqe-task");

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
	for (auto r : result) {
		std::stringstream ss;
		ss << std::setprecision(12) << r.second << " at (" << r.first.transpose() << ")";
		XACCInfo("Computed VQE Energy = " + ss.str());
	}

	xacc::Finalize();
}




