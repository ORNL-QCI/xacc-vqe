#include "XACC.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include "PauliOperator.hpp"
#include "VQEProgram.hpp"
#include "VQETask.hpp"
#include "VQEParameterGenerator.hpp"
#include "DiagonalizeTask.hpp"
#include "FermionToSpinTransformation.hpp"
#include "GateFunction.hpp"

namespace py = pybind11;

using namespace xacc;
using namespace xacc::vqe;

using GateFunctionPtr = std::shared_ptr<xacc::quantum::GateFunction>;

PauliOperator compileFermionHamiltonian(const std::string& fermiSrc) {

	auto mpi4py = pybind11::module::import("mpi4py.MPI");
	auto comm = mpi4py.attr("COMM_WORLD");

	if (!xacc::isInitialized()) {
		xacc::Initialize();
		XACCInfo("You did not initialize the XACC framework. "
				"Auto-running xacc::Initialize().");
	}

	// Get the user-specified Accelerator,
	// or TNQVM if none specified
	xacc::setAccelerator("vqe-dummy");
	// Set the default Accelerator to TNQVM
	if (xacc::hasAccelerator("tnqvm")) {
		xacc::setAccelerator("tnqvm");
	}

	auto accelerator = xacc::getAccelerator();
	boost::mpi::communicator world;

	xacc::setOption("vqe-task", "vqe-profile");

	auto program = std::make_shared<VQEProgram>(accelerator, fermiSrc, world);
	program->build();

	auto parameters = VQEParameterGenerator::generateParameters(program->getNParameters(), world);

	return program->getPauliOperator();
}

/**
 * kwargs can be {'accelerator':'name', 'task':'name', 'ansatz':GateFunctionPtr,
 * 'diagonalize-backend':'backendname', 'n-qubits':nqubitsint }
 * @param op
 * @param kwargs
 * @return
 */
VQETaskResult execute(PauliOperator& op, py::kwargs kwargs) {

	auto mpi4py = pybind11::module::import("mpi4py.MPI");
	auto comm = mpi4py.attr("COMM_WORLD");

	if (!xacc::isInitialized()) {
		xacc::Initialize();
		XACCInfo("You did not initialize the XACC framework. "
				"Auto-running xacc::Initialize().");
	}

	// Get the user-specified Accelerator,
	// or TNQVM if none specified
	xacc::setAccelerator("vqe-dummy");
	// Set the default Accelerator to TNQVM
	if (xacc::hasAccelerator("tnqvm")) {
		xacc::setAccelerator("tnqvm");
	}

	auto accelerator = xacc::getAccelerator();
	auto statePrep = GateFunctionPtr(nullptr);
	boost::mpi::communicator world;

	// Get the task to run
	std::string task = "vqe-diagonalize";

	int maxQbit = 0;
	for (auto& kv : op.getTerms()) {
		if (!kv.second.isIdentity()) {
			int qbit = kv.second.ops().rbegin()->first;
			if (qbit > maxQbit)
				maxQbit = qbit;
		}
	}
	int nQubits = maxQbit + 1;

	if (kwargs) {
		task = kwargs.contains("task") ?
				kwargs["task"].cast<std::string>() : task;
		std::string accStr =
				kwargs.contains("accelerator") ?
						kwargs["accelerator"].cast<std::string>() : "";
		if (!accStr.empty())
			accelerator = xacc::getAccelerator(accStr);

		if (kwargs.contains("diagonalize-backend")) {
			xacc::setOption("diagonalize-backend",
					kwargs["diagonalize-backend"].cast<std::string>());
		}

		nQubits =
				kwargs.contains("n-qubits") ?
						kwargs["n-qubits"].cast<int>() : nQubits;

		if (kwargs.contains("ansatz")) {
			statePrep = kwargs["ansatz"].cast<GateFunctionPtr>();
		}

		if (kwargs.contains("vqe-params")) {
			xacc::setOption("vqe-parameters", kwargs["vqe-params"].cast<std::string>());
		}
	}

	XACCInfo("XACC VQE Python set n-qubits = " + std::to_string(nQubits));

	xacc::setOption("vqe-task", task);
	xacc::setOption("n-qubits", std::to_string(nQubits));

	auto program = std::make_shared<VQEProgram>(accelerator, op, statePrep, world);
	program->build();

	auto parameters = VQEParameterGenerator::generateParameters(program->getNParameters(), world);
	auto vqeTask = xacc::ServiceRegistry::instance()->getService<VQETask>(task);
	vqeTask->setVQEProgram(program);

	return vqeTask->execute(parameters);
}

VQETaskResult execute(py::object& fermionOperator, py::kwargs kwargs) {

	if (!fermionOperator) {
		XACCError("FermionOperator was null. Exiting.");
	}

	auto mpi4py = pybind11::module::import("mpi4py.MPI");
	auto comm = mpi4py.attr("COMM_WORLD");

	if(!pybind11::hasattr(fermionOperator, "terms")) {
		XACCError("This is not a FermionOperator, it does not have a terms dict");
	}

	if (!xacc::isInitialized()) {
		xacc::Initialize();
		XACCInfo("You did not initialize the XACC framework. "
				"Auto-running xacc::Initialize().");
	}

	auto terms = fermionOperator.attr("terms").cast<py::dict>();

	std::stringstream s;
	s << "__qpu__ openfermion_kernel() {\n";
	for (auto& kv : terms) {

		auto fTerm = kv.first.cast<py::tuple>();
		auto coeff = kv.second.cast<std::complex<double>>();
		s << std::real(coeff) << " ";
		for (auto t : fTerm) {
			auto ct = t.cast<py::tuple>();
			s << ct[0].cast<int>() << " " << ct[1].cast<int>() << " ";
		}
		s << "\n";
	}
	s << "}";

	// Get the user-specified Accelerator,
	// or TNQVM if none specified
	xacc::setAccelerator("vqe-dummy");
	// Set the default Accelerator to TNQVM
	if (xacc::hasAccelerator("tnqvm")) {
		xacc::setAccelerator("tnqvm");
	}

	auto accelerator = xacc::getAccelerator();
	auto statePrep = GateFunctionPtr(nullptr);
	boost::mpi::communicator world;

	// Get the task to run
	std::string task = "vqe-diagonalize";

	if (kwargs) {
		task = kwargs.contains("task") ?
				kwargs["task"].cast<std::string>() : task;
		std::string accStr =
				kwargs.contains("accelerator") ?
						kwargs["accelerator"].cast<std::string>() : "";
		if (!accStr.empty())
			accelerator = xacc::getAccelerator(accStr);

		if (kwargs.contains("diagonalize-backend")) {
			xacc::setOption("diagonalize-backend",
					kwargs["diagonalize-backend"].cast<std::string>());
		}

		if (kwargs.contains("ansatz")) {
			statePrep = kwargs["ansatz"].cast<GateFunctionPtr>();
		}

		if (kwargs.contains("vqe-params")) {
			xacc::setOption("vqe-parameters",
					kwargs["vqe-params"].cast<std::string>());
		}
	}

	xacc::setOption("vqe-task", task);

	auto program = std::make_shared<VQEProgram>(accelerator, s.str(), statePrep,
			world);
	program->build();

	auto parameters = VQEParameterGenerator::generateParameters(
			program->getNParameters(), world);
	auto vqeTask = xacc::ServiceRegistry::instance()->getService<VQETask>(task);
	vqeTask->setVQEProgram(program);

	return vqeTask->execute(parameters);
}

PYBIND11_MODULE(pyxaccvqe, m) {
	m.doc() = "Python bindings for XACC VQE.";

	py::add_ostream_redirect(m, "ostream_redirect");

	py::class_<VQETaskResult>(m, "VQETaskResult")
	    .def_readonly("buffers", &VQETaskResult::buffers)
		.def_readonly("results", &VQETaskResult::results);

	py::class_<PauliOperator>(m,"PauliOperator")
			.def(py::init<>())
			.def(py::init<std::complex<double>>())
			.def(py::init<double>())
			.def(py::init<std::string>())
			.def(py::init<std::map<int, std::string>>())
			.def(py::init<std::map<int, std::string>, double>())
			.def(py::init<std::map<int, std::string>, std::complex<double>>())
			.def(py::init<std::map<int, std::string>, std::string>())
			.def(py::init<std::map<int, std::string>, std::complex<double>, std::string>())
			.def(py::self + py::self)
			.def(py::self += py::self)
			.def(py::self *= py::self)
			.def(py::self *= double())
			.def(py::self * py::self)
			.def(py::self *= std::complex<double>())
			.def(py::self -= py::self)
			.def(py::self - py::self)
			.def("__repr__", &PauliOperator::toString)
			.def("eval", &PauliOperator::eval)
			.def("toXACCIR", &PauliOperator::toXACCIR);


	m.def("execute",
			[](PauliOperator& op, py::kwargs kwargs) -> VQETaskResult {
				py::scoped_ostream_redirect stream(
						std::cout,                              // std::ostream&
						py::module::import("sys").attr("stdout")// Python output
				);
				return execute(op, kwargs);
			});

	m.def("execute",
			[](py::object& op, py::kwargs kwargs) -> VQETaskResult {
				py::scoped_ostream_redirect stream(
						std::cout,                              // std::ostream&
						py::module::import("sys").attr("stdout")// Python output
				);
				return execute(op, kwargs);
			});

	m.def("compileFermionHamiltonian",
			[](const std::string& src) -> PauliOperator {
				py::scoped_ostream_redirect stream(
						std::cout,                              // std::ostream&
						py::module::import("sys").attr("stdout")// Python output
				);
				return compileFermionHamiltonian(src);
			});
}

