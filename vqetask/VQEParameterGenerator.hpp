#ifndef VQETASK_VQEPARAMETERGENERATOR_HPP_
#define VQETASK_VQEPARAMETERGENERATOR_HPP_

#include "XACC.hpp"
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/mpi.hpp>

namespace xacc {
namespace vqe {

class VQEParameterGenerator {

public:

	static Eigen::VectorXd generateParameters(const int nParameters, boost::mpi::communicator& comm) {

		if (xacc::optionExists("vqe-parameters")) {
			auto paramStr = xacc::getOption("vqe-parameters");
			std::vector<std::string> split;
			boost::split(split, paramStr, boost::is_any_of(","));

			if (xacc::getOption("vqe-task") == "sweep-1d") {
				return Eigen::VectorXd::LinSpaced(50, std::stod(split[0]), std::stod(split[1]));
			} else {
				Eigen::VectorXd params(nParameters);
				for (int i = 0; i < split.size(); i++) {
					params(i) = std::stod(split[i]);
				}

				return params;
			}
		} else {
			std::srand(time(0));
			auto pi = boost::math::constants::pi<double>();
			Eigen::VectorXd rand;

			std::vector<double> data;

			// Random parameters between -pi and pi
			if (comm.rank() == 0) {
				rand = -1.0 * pi * Eigen::VectorXd::Ones(nParameters)
						+ (Eigen::VectorXd::Random(nParameters) * 0.5
								+ Eigen::VectorXd::Ones(nParameters) * 0.5)
								* (pi - (-1 * pi));
				data.resize(rand.size());
				Eigen::VectorXd::Map(&data[0], rand.size()) = rand;
			}

			boost::mpi::broadcast(comm, data, 0);

			if (comm.rank() != 0) {
				rand = Eigen::Map<Eigen::VectorXd>(data.data(), data.size());
			}
			return rand;
		}
	}

};
}
}


#endif /* VQETASK_STATEPREPARATIONEVALUATOR_HPP_ */
