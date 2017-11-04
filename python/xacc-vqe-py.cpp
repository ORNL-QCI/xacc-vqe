#include "XACC.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>

#include "VQEProblem.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyxaccvqe, m) {
    m.doc() = "Python bindings for XACC VQE.";

	// Expose the Program object
	py::class_<xacc::vqe::VQEProblem> problem(m,
			"VQEProblem", "");
	problem.def(py::init<const std::string &>(), "The constructor");
	problem.def(py::init<const std::string &, const int>(), "The constructor");

	problem.def("minimize", (const Eigen::VectorXd (xacc::vqe::VQEProblem::*)()) &xacc::vqe::VQEProblem::minimize, "");

	problem.def("minimize", (const Eigen::VectorXd (xacc::vqe::VQEProblem::*)(Eigen::VectorXd&)) &xacc::vqe::VQEProblem::minimize,
			"");
	problem.def("getCurrentEnergy", &xacc::vqe::VQEProblem::getCurrentEnergy, "Run VQE.");

}


