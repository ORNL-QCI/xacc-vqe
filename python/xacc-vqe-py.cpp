#include "XACC.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include "VQEProblem.hpp"

namespace py = pybind11;

std::pair<Eigen::VectorXd, double> execute_vqe(const std::string& hamiltonianSrc,
		const std::string& statePrepSrc, const int nQbits, const std::string& accelerator = "tnqvm") {
	std::ofstream tempStatePrep(".tempStatePrep.hpp");
	tempStatePrep << statePrepSrc;
	tempStatePrep.close();
	xacc::setOption("vqe-state-prep-kernel", ".tempStatePrep.hpp");
	xacc::setOption("n-qubits", std::to_string(nQbits));
	xacc::setOption("accelerator", accelerator);

	xacc::vqe::VQEProblem problem(hamiltonianSrc);
	auto params = problem.minimize();
	std::remove(".tempStatePrep.hpp");

	return std::make_pair(params, problem.currentEnergy);
}

std::pair<Eigen::VectorXd, double> execute_vqe(const std::string& hamiltonianSrc,
                const std::string& statePrepSrc, const int nQbits, Eigen::VectorXd& params, const std::string& accelerator = "tnqvm") {
        std::ofstream tempStatePrep(".tempStatePrep.hpp");
        tempStatePrep << statePrepSrc;
        tempStatePrep.close();
        xacc::setOption("vqe-state-prep-kernel", ".tempStatePrep.hpp");
        xacc::setOption("n-qubits", std::to_string(nQbits));
        xacc::setOption("accelerator", accelerator);

        xacc::vqe::VQEProblem problem(hamiltonianSrc);
        auto newParams = problem.minimize(params);
        std::remove(".tempStatePrep.hpp");

        return std::make_pair(newParams, problem.currentEnergy);
}

std::vector<double> sweepParameter1D(const std::string& hamiltonianSrc,
                const std::string& statePrepSrc, const int nQbits, const Eigen::VectorXd& range, const std::string& accelerator = "tnqvm") {
    std::ofstream tempStatePrep(".tempStatePrep.hpp");
    tempStatePrep << statePrepSrc;
    tempStatePrep.close();
    xacc::setOption("vqe-state-prep-kernel", ".tempStatePrep.hpp");
    xacc::setOption("n-qubits", std::to_string(nQbits));
    xacc::setOption("accelerator", accelerator);
    xacc::vqe::VQEProblem problem(hamiltonianSrc);

    std::vector<double> energies;
    Eigen::VectorXd params(1);
    for (int i = 0; i < range.rows(); i++) {
    	params(0) = range(i);
    	energies.push_back(problem(params));
    }
    return energies;
}

double computeEnergyAtParameters(
		const std::string& hamiltonianSrc, const std::string& statePrepSrc,
		const int nQbits, Eigen::VectorXd& params,
		const std::string& accelerator = "tnqvm") {

	std::ofstream tempStatePrep(".tempStatePrep.hpp");
    tempStatePrep << statePrepSrc;
    tempStatePrep.close();
    xacc::setOption("vqe-state-prep-kernel", ".tempStatePrep.hpp");
    xacc::setOption("n-qubits", std::to_string(nQbits));
    xacc::setOption("accelerator", accelerator);
    xacc::vqe::VQEProblem problem(hamiltonianSrc);

    auto energy = problem(params);

    std::cout << "I HAVE THE ENERGY AT " << energy << "\n";

    std::remove(".tempStatePrep.hpp");

    return energy;
}

PYBIND11_MODULE(pyxaccvqe, m) {
    m.doc() = "Python bindings for XACC VQE.";

py::add_ostream_redirect(m, "ostream_redirect");

    m.def("execute_vqe", [](const std::string& hamiltonianSrc,
    		const std::string& statePrepSrc, const int nQbits, const std::string& accelerator = "tnqvm") -> std::pair<Eigen::VectorXd, double> {
        py::scoped_ostream_redirect stream(
            std::cout,                               // std::ostream&
            py::module::import("sys").attr("stdout") // Python output
        );
        return execute_vqe(hamiltonianSrc, statePrepSrc, nQbits, accelerator);
    });

    m.def("execute_vqe", [](const std::string& hamiltonianSrc,
                const std::string& statePrepSrc, const int nQbits, Eigen::VectorXd& params, const std::string& accelerator = "tnqvm") -> std::pair<Eigen::VectorXd, double> {
        py::scoped_ostream_redirect stream(
            std::cout,                               // std::ostream&
            py::module::import("sys").attr("stdout") // Python output
        );
        return execute_vqe(hamiltonianSrc, statePrepSrc, nQbits, params, accelerator);
    });

    m.def("sweepParameter1D", [](const std::string& hamiltonianSrc,
            const std::string& statePrepSrc, const int nQbits, const Eigen::VectorXd& range, const std::string& accelerator = "tnqvm") -> std::vector<double> {
          py::scoped_ostream_redirect stream(
              std::cout,                               // std::ostream&
              py::module::import("sys").attr("stdout") // Python output
          );
          return sweepParameter1D(hamiltonianSrc, statePrepSrc, nQbits, range, accelerator);
      });

 m.def("computeEnergyAtParameters", [](
    		const std::string& hamiltonianSrc, const std::string& statePrepSrc,
    		const int nQbits, Eigen::VectorXd& params,
    		const std::string& accelerator = "tnqvm") -> double {
          py::scoped_ostream_redirect stream(
              std::cout,                               // std::ostream&
              py::module::import("sys").attr("stdout") // Python output
          );
          return computeEnergyAtParameters(hamiltonianSrc, statePrepSrc, nQbits, params, accelerator);
      });

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


