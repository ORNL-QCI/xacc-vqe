#include "XACC.hpp"
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>
// #include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "AcceleratorDecorator.hpp"

#include "MPIProvider.hpp"

#include "GateFunction.hpp"
#include "PauliOperator.hpp"
#include "VQEParameterGenerator.hpp"
#include "VQEProgram.hpp"
#include "VQETask.hpp"

namespace py = pybind11;

using namespace xacc;
using namespace xacc::vqe;

using GateFunctionPtr = std::shared_ptr<xacc::Function>;

/**
 * Compile the given source code string and produce
 * the corresponding PauliOperator instance.
 *
 * @param fermiSrc Source code describing Fermionic Hamiltonian
 * @return op Pauli Hamiltonian representation of given Fermion Hamiltonian
 * source code
 *
 */

PauliOperator compile(const std::string &fermiSrc) {

  if (!xacc::isInitialized()) {
    xacc::Initialize({"--use-cout", "--no-color"});
    xacc::info("You did not initialize the XACC framework. "
               "Auto-running xacc::Initialize().");
  }

  std::shared_ptr<MPIProvider> provider;
  if (xacc::hasService<MPIProvider>("boost-mpi")) {
    provider = xacc::getService<MPIProvider>("boost-mpi");
    auto mpi4py = pybind11::module::import("mpi4py.MPI");
    auto comm = mpi4py.attr("COMM_WORLD");
  } else {
    provider = xacc::getService<MPIProvider>("no-mpi");
  }

  provider->initialize();
  auto world = provider->getCommunicator();

  // Get the user-specified Accelerator,
  // or TNQVM if none specified
  xacc::setAccelerator("vqe-dummy");
  // Set the default Accelerator to TNQVM
  if (xacc::hasAccelerator("tnqvm")) {
    xacc::setAccelerator("tnqvm");
  }

  auto accelerator = xacc::getAccelerator();

  // Set to vqe-profile because it doesn't require state prep
  xacc::setOption("vqe-task", "vqe-profile");
  auto program = std::make_shared<VQEProgram>(accelerator, fermiSrc, world);
  program->build();

  //	xacc::clearOptions();
  return program->getPauliOperator();
}

PauliOperator compile(py::object fermionOperator, py::kwargs kwargs) {

  if (!xacc::isInitialized()) {
    xacc::Initialize({"--use-cout", "--no-color"});
    xacc::info("You did not initialize the XACC framework. "
               "Auto-running xacc::Initialize().");
  }

  std::shared_ptr<MPIProvider> provider;
  if (xacc::hasService<MPIProvider>("boost-mpi")) {
    provider = xacc::getService<MPIProvider>("boost-mpi");
    auto mpi4py = pybind11::module::import("mpi4py.MPI");
    auto comm = mpi4py.attr("COMM_WORLD");
  } else {
    provider = xacc::getService<MPIProvider>("no-mpi");
  }

  provider->initialize();
  auto world = provider->getCommunicator();

  auto terms = fermionOperator.attr("terms").cast<py::dict>();

  std::stringstream s;
  s << "__qpu__ openfermion_kernel() {\n";
  for (auto &kv : terms) {

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
  if (!xacc::optionExists("accelerator")) {
    xacc::setAccelerator("vqe-dummy");
    // Set the default Accelerator to TNQVM
    if (xacc::hasAccelerator("tnqvm")) {
      xacc::setAccelerator("tnqvm");
    }
  }

  auto accelerator = xacc::getAccelerator();

  // Set to vqe-profile because it doesn't require state prep
  xacc::setOption("vqe-task", "vqe-profile");
  auto program = std::make_shared<VQEProgram>(accelerator, s.str(), world);
  program->build();
  //	xacc::clearOptions();
  return program->getPauliOperator();
}

/**
 * kwargs can be {'accelerator':'name', 'task':'name', 'ansatz':GateFunctionPtr,
 * 'diagonalize-backend':'backendname', 'n-qubits':nqubitsint }
 * @param op
 * @param kwargs
 * @return
 */
VQETaskResult execute(PauliOperator &op,
                      std::shared_ptr<AcceleratorBuffer> buffer,
                      py::kwargs kwargs) {
  if (!xacc::isInitialized()) {
    xacc::Initialize({"--use-cout", "--no-color"});
    xacc::info("You did not initialize the XACC framework. "
               "Auto-running xacc::Initialize().");
  }

  std::shared_ptr<MPIProvider> provider;
  if (kwargs.contains("mpi-provider")) {
    provider = xacc::getService<MPIProvider>(
        kwargs["mpi-provider"].cast<std::string>());
  } else {
    if (xacc::hasService<MPIProvider>("boost-mpi")) {
      provider = xacc::getService<MPIProvider>("boost-mpi");
      auto mpi4py = pybind11::module::import("mpi4py.MPI");
      auto comm = mpi4py.attr("COMM_WORLD");
    } else {
      provider = xacc::getService<MPIProvider>("no-mpi");
    }
  }

  provider->initialize();
  auto world = provider->getCommunicator();

  // Get the user-specified Accelerator,
  // or TNQVM if none specified
  if (!xacc::optionExists("accelerator")) {
    xacc::setAccelerator("vqe-dummy");
    // Set the default Accelerator to TNQVM
    if (xacc::hasAccelerator("tnqvm")) {
      xacc::setAccelerator("tnqvm");
    }
  }

  auto accelerator = xacc::getAccelerator();
  auto statePrep = GateFunctionPtr(nullptr);

  // Get the task to run
  std::string task = "vqe-diagonalize";

  int maxQbit = 0;
  for (auto &kv : op.getTerms()) {
    if (!kv.second.isIdentity()) {
      int qbit = kv.second.ops().rbegin()->first;
      if (qbit > maxQbit)
        maxQbit = qbit;
    }
  }
  int nQubits = maxQbit + 1;

  if (kwargs) {
    task = kwargs.contains("task") ? kwargs["task"].cast<std::string>() : task;

    if (kwargs.contains("accelerator")) {
      auto obj = kwargs["accelerator"];
      if (py::isinstance<std::shared_ptr<Accelerator>>(obj) ||
          py::isinstance<std::shared_ptr<AcceleratorDecorator>>(obj) ||
          py::isinstance<AcceleratorDecorator>(obj) ||
          py::isinstance<Accelerator>(obj)) {
        accelerator = obj.cast<std::shared_ptr<Accelerator>>();
      } else {
        std::string accStr = kwargs.contains("accelerator")
                                 ? kwargs["accelerator"].cast<std::string>()
                                 : "";
        if (!accStr.empty())
          accelerator = xacc::getAccelerator(accStr);
      }
    }

    if (kwargs.contains("diagonalize-backend")) {
      xacc::setOption("diagonalize-backend",
                      kwargs["diagonalize-backend"].cast<std::string>());
    }

    nQubits =
        kwargs.contains("n-qubits") ? kwargs["n-qubits"].cast<int>() : nQubits;

    if (kwargs.contains("ansatz")) {
      statePrep = kwargs["ansatz"].cast<std::shared_ptr<Function>>();
    }

    if (kwargs.contains("vqe-params")) {
      xacc::setOption("vqe-parameters",
                      kwargs["vqe-params"].cast<std::string>());
    }

    if (kwargs.contains("n-electrons")) {
      auto nElectrons = kwargs["n-electrons"].cast<int>();
      xacc::setOption("n-electrons", std::to_string(nElectrons));
    }

    if (kwargs.contains("error-mitigation")) {
      auto errorMitigationStrategies =
          kwargs["error-mitigation"].cast<std::vector<std::string>>();
      for (auto e : errorMitigationStrategies) {
        xacc::setOption(e, "");
      }
    }

    if (kwargs.contains("qubit-map")) {
      auto qbitmap = kwargs["qubit-map"].cast<std::vector<int>>();
      auto mapStr = std::to_string(qbitmap[0]);
      for (int i = 1; i < qbitmap.size(); i++) {
        mapStr += "," + std::to_string(qbitmap[i]);
      }
      xacc::setOption("qubit-map", mapStr);
    }
  }

  xacc::setOption("vqe-task", task);
  xacc::setOption("n-qubits", std::to_string(nQubits));

  auto program =
      std::make_shared<VQEProgram>(accelerator, op, statePrep, world);
  program->build();
  program->setGlobalBuffer(buffer);

  auto parameters = VQEParameterGenerator::generateParameters(
      program->getNParameters(), world);
  auto vqeTask = xacc::getService<VQETask>(task);
  vqeTask->setVQEProgram(program);
  auto result = vqeTask->execute(parameters);
  //	xacc::clearOptions();
  return result;
}

VQETaskResult execute(py::object &fermionOperator,
                      std::shared_ptr<AcceleratorBuffer> buffer,
                      py::kwargs kwargs) {

  if (!fermionOperator) {
    xacc::error("FermionOperator was null. Exiting.");
  }

  if (!pybind11::hasattr(fermionOperator, "terms")) {
    xacc::error("This is not a FermionOperator, it does not have a terms dict");
  }

  if (!xacc::isInitialized()) {
    xacc::Initialize({"--use-cout", "--no-color"});
    xacc::info("You did not initialize the XACC framework. "
               "Auto-running xacc::Initialize().");
  }

  std::shared_ptr<MPIProvider> provider;
  if (xacc::hasService<MPIProvider>("boost-mpi")) {
    provider = xacc::getService<MPIProvider>("boost-mpi");
    auto mpi4py = pybind11::module::import("mpi4py.MPI");
    auto comm = mpi4py.attr("COMM_WORLD");
  } else {
    provider = xacc::getService<MPIProvider>("no-mpi");
  }

  provider->initialize();
  auto world = provider->getCommunicator();

  auto terms = fermionOperator.attr("terms").cast<py::dict>();

  std::stringstream s;
  s << "__qpu__ openfermion_kernel() {\n";
  for (auto &kv : terms) {

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
  if (!xacc::optionExists("accelerator")) {
    xacc::setAccelerator("vqe-dummy");
    // Set the default Accelerator to TNQVM
    if (xacc::hasAccelerator("tnqvm")) {
      xacc::setAccelerator("tnqvm");
    }
  }

  auto accelerator = xacc::getAccelerator();
  auto statePrep = GateFunctionPtr(nullptr);

  // Get the task to run
  std::string task = "vqe-diagonalize";

  if (kwargs) {
    task = kwargs.contains("task") ? kwargs["task"].cast<std::string>() : task;
    if (kwargs.contains("accelerator")) {
      auto obj = kwargs["accelerator"];
      if (py::isinstance<std::shared_ptr<Accelerator>>(obj) ||
          py::isinstance<std::shared_ptr<AcceleratorDecorator>>(obj) ||
          py::isinstance<AcceleratorDecorator>(obj) ||
          py::isinstance<Accelerator>(obj)) {
        accelerator = obj.cast<std::shared_ptr<Accelerator>>();
      } else {
        std::string accStr = kwargs.contains("accelerator")
                                 ? kwargs["accelerator"].cast<std::string>()
                                 : "";
        if (!accStr.empty())
          accelerator = xacc::getAccelerator(accStr);
      }
    }

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

    if (kwargs.contains("n-electrons")) {
      auto nElectrons = kwargs["n-electrons"].cast<int>();
      xacc::setOption("n-electrons", std::to_string(nElectrons));
    }

    if (kwargs.contains("error-mitigation")) {
      auto errorMitigationStrategies =
          kwargs["error-mitigation"].cast<std::vector<std::string>>();
      for (auto e : errorMitigationStrategies) {
        xacc::setOption(e, "");
      }
    }

    if (kwargs.contains("qubit-map")) {
      auto qbitmap = kwargs["qubit-map"].cast<std::vector<int>>();
      auto mapStr = std::to_string(qbitmap[0]);
      for (int i = 1; i < qbitmap.size(); i++) {
        mapStr += "," + std::to_string(qbitmap[i]);
      }
      xacc::setOption("qubit-map", mapStr);
    }

    if (kwargs.contains("transformation")) {
      xacc::setOption("fermion-transformation",
                      kwargs["transformation"].cast<std::string>());
    }
  }

  xacc::setOption("vqe-task", task);

  auto program =
      std::make_shared<VQEProgram>(accelerator, s.str(), statePrep, world);
  program->build();
  program->setGlobalBuffer(buffer);

  auto parameters = VQEParameterGenerator::generateParameters(
      program->getNParameters(), world);
  auto vqeTask = xacc::getService<VQETask>(task);
  vqeTask->setVQEProgram(program);

  auto result = vqeTask->execute(parameters);

  //	xacc::clearOptions();
  return result;
}

PYBIND11_MODULE(_pyxaccvqe, m) {
  m.doc() = "Python bindings for XACC VQE.";

  py::add_ostream_redirect(m, "ostream_redirect");

  py::class_<VQETaskResult>(m, "VQETaskResult")
      .def(py::init<double, Eigen::VectorXd>())
      .def_readonly("angles", &VQETaskResult::angles)
      .def_readonly("nQpuCalls", &VQETaskResult::nQpuCalls)
      .def_readonly("vqeIterations", &VQETaskResult::vqeIterations)
      .def_readonly("energy", &VQETaskResult::energy)
      .def_readonly("ansatzQASM", &VQETaskResult::ansatzQASM);

  py::class_<Term>(m, "Term").def("coeff", &Term::coeff).def("ops",&Term::ops);
  py::class_<PauliOperator>(m, "PauliOperator")
      .def(py::init<>())
      .def(py::init<std::complex<double>>())
      .def(py::init<double>())
      .def(py::init<std::string>())
      .def(py::init<std::map<int, std::string>>())
      .def(py::init<std::map<int, std::string>, double>())
      .def(py::init<std::map<int, std::string>, std::complex<double>>())
      .def(py::init<std::map<int, std::string>, std::string>())
      .def(py::init<std::map<int, std::string>, std::complex<double>,
                    std::string>())
      .def(py::self + py::self)
      .def(py::self += py::self)
      .def(py::self *= py::self)
      .def(py::self *= double())
      .def(py::self * py::self)
      .def(py::self *= std::complex<double>())
      .def(py::self -= py::self)
      .def(py::self - py::self)
      .def("__eq__", &PauliOperator::operator==)
      .def("__repr__", &PauliOperator::toString)
      .def("eval", &PauliOperator::eval)
      .def("toXACCIR", &PauliOperator::toXACCIR)
      .def("fromXACCIR", &PauliOperator::fromXACCIR)
      .def("fromString", &PauliOperator::fromString)
      .def("nTerms", &PauliOperator::nTerms)
      .def("isClose", &PauliOperator::isClose)
      .def("commutes", &PauliOperator::commutes)
      .def("__len__", &PauliOperator::nTerms)
      .def("nQubits", &PauliOperator::nQubits)
      .def("__iter__",
           [](PauliOperator &op) {
             return py::make_iterator(op.begin(), op.end());
           },
           py::keep_alive<0, 1>());

  m.def("execute",
        (VQETaskResult(*)(PauliOperator & op,
                          std::shared_ptr<AcceleratorBuffer> b,
                          py::kwargs kwargs)) &
            execute,
        py::call_guard<py::scoped_ostream_redirect,
                       py::scoped_estream_redirect>(),
        "");

  m.def("execute",
        (VQETaskResult(*)(py::object & op, std::shared_ptr<AcceleratorBuffer> b,
                          py::kwargs kwargs)) &
            execute,
        py::call_guard<py::scoped_ostream_redirect,
                       py::scoped_estream_redirect>(),
        "");

  m.def("compile", (PauliOperator(*)(const std::string &src)) & compile,
        py::call_guard<py::scoped_ostream_redirect,
                       py::scoped_estream_redirect>(),
        "");

  m.def("compile", [](py::object &op, py::kwargs kwargs) -> PauliOperator {
    py::scoped_ostream_redirect stream(
        std::cout,                               // std::ostream&
        py::module::import("sys").attr("stdout") // Python output
    );
    return compile(op, kwargs);
  });
}
