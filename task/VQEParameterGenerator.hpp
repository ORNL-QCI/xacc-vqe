#ifndef TASK_VQEPARAMETERGENERATOR_HPP_
#define TASK_VQEPARAMETERGENERATOR_HPP_

#include "XACC.hpp"
#include <Eigen/Dense>

namespace xacc {
namespace vqe {

class VQEParameterGenerator {

public:

	static Eigen::VectorXd generateParameters(const int nParameters, std::shared_ptr<Communicator> comm) {

		if (xacc::optionExists("vqe-parameters")) {
			if (xacc::getOption("vqe-task") == "sweep-1d") {
				auto paramStr = xacc::getOption("vqe-parameters");

				// HERE WE COULD HAVE SOMETHING LIKE EITHER
				// 50:-3.14,3.14 or
				// -3.14,3.14
				std::vector<std::string> split, colon;
                split = xacc::split(paramStr, ',');
				// boost::split(split, paramStr, boost::is_any_of(","));
				int nSteps = 50;
				std::string minVal;
				if (paramStr.find(":") != std::string::npos) {
					colon = xacc::split(split[0], ':');
					nSteps = std::stoi(colon[0]);
					minVal = colon[1];
				} else {
					minVal = split[0];
				}

				return Eigen::VectorXd::LinSpaced(nSteps, std::stod(minVal),
						std::stod(split[1]));
			} else {
				auto paramStr = xacc::getOption("vqe-parameters");
				std::vector<std::string> split;
				split = xacc::split( paramStr, ',');
				Eigen::VectorXd params(nParameters);
				for (int i = 0; i < split.size(); i++) {
					params(i) = std::stod(split[i]);
				}

				return params;
			}
		} else {
			std::srand(time(0));
			auto pi = 3.141592653589793238;
			Eigen::VectorXd rand;

			std::vector<double> data;

			// Random parameters between -pi and pi
			if (comm->rank() == 0) {
				rand = -1.0 * pi * Eigen::VectorXd::Ones(nParameters)
						+ (Eigen::VectorXd::Random(nParameters) * 0.5
								+ Eigen::VectorXd::Ones(nParameters) * 0.5)
								* (pi - (-1 * pi));
				data.resize(rand.size());
				Eigen::VectorXd::Map(&data[0], rand.size()) = rand;
			}

			comm->broadcast(data, 0);

			if (comm->rank() != 0) {
				rand = Eigen::Map<Eigen::VectorXd>(data.data(), data.size());
			}
			return rand;
		}
	}

};
}
}


#endif /* VQETASK_STATEPREPARATIONEVALUATOR_HPP_ */
