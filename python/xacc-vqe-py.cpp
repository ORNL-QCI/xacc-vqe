#include "XACC.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include "VQEProblem.hpp"
#include "VQEProgram.hpp"
#include "VQETask.hpp"
#include "VQEParameterGenerator.hpp"

namespace py = pybind11;

using namespace xacc;
using namespace xacc::vqe;

std::pair<Eigen::VectorXd, double> execute_vqe(
		const std::string& hamiltonianSrc, const std::string& statePrepSrc,
		const int nQbits, const std::string& accelerator = "tnqvm") {

	xacc::setOption("n-qubits", std::to_string(nQbits));
	auto qpu = xacc::getAccelerator(accelerator);

	boost::mpi::communicator comm;
	auto program = std::make_shared<VQEProgram>(qpu, hamiltonianSrc, statePrepSrc, comm);
	program->build();

	auto task = ServiceRegistry::instance()->getService<VQETask>("vqe");
	task->setVQEProgram(program);

	auto parameters = VQEParameterGenerator::generateParameters(program->getNParameters(), comm);
	auto result = task->execute(parameters);
	return result[0];
}

std::pair<Eigen::VectorXd, double> execute_vqe(
		const std::string& hamiltonianSrc, const std::string& statePrepSrc,
		const int nQbits, Eigen::VectorXd& params,
		const std::string& accelerator = "tnqvm") {
	xacc::setOption("n-qubits", std::to_string(nQbits));
	auto qpu = xacc::getAccelerator(accelerator);

	boost::mpi::communicator comm;
	auto program = std::make_shared<VQEProgram>(qpu, hamiltonianSrc, statePrepSrc, comm);
	program->build();

	auto task = ServiceRegistry::instance()->getService<VQETask>("vqe");
	task->setVQEProgram(program);

	auto result = task->execute(params);
	return result[0];
}

std::vector<double> sweepParameter1D(const std::string& hamiltonianSrc,
		const std::string& statePrepSrc, const int nQbits,
		const Eigen::VectorXd& range,
		const std::string& accelerator = "tnqvm") {

	xacc::setOption("n-qubits", std::to_string(nQbits));
	auto qpu = xacc::getAccelerator(accelerator);

	boost::mpi::communicator comm;
	auto program = std::make_shared<VQEProgram>(qpu, hamiltonianSrc, statePrepSrc, comm);
	program->build();

	auto task = ServiceRegistry::instance()->getService<VQETask>("sweep-1d");
	task->setVQEProgram(program);

	auto result = task->execute(range);

	std::vector<double> energies;
	for (auto r : result) {
		energies.push_back(r.second);
	}

    return energies;
}

double computeEnergyAtParameters(
		const std::string& hamiltonianSrc, const std::string& statePrepSrc,
		const int nQbits, Eigen::VectorXd& params,
		const std::string& accelerator = "tnqvm") {

	xacc::setOption("n-qubits", std::to_string(nQbits));
	auto qpu = xacc::getAccelerator(accelerator);

	boost::mpi::communicator comm;
	auto program = std::make_shared<VQEProgram>(qpu, hamiltonianSrc,
			statePrepSrc, comm);
	program->build();

	auto task = ServiceRegistry::instance()->getService<VQETask>(
			"compute-energy");
	task->setVQEProgram(program);

	auto result = task->execute(params);
	return result[0].second;
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
}


