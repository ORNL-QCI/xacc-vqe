/*
 * CommutingSetGenerator.hpp
 *
 *  Created on: Aug 4, 2017
 *      Author: aqw
 */

#ifndef VQE_TRANSFORMATION_COMMUTINGSETGENERATOR_HPP_
#define VQE_TRANSFORMATION_COMMUTINGSETGENERATOR_HPP_

#include "SpinInstruction.hpp"
#include <Eigen/Core>
#include <numeric>

namespace xacc {

namespace vqe {
class CommutingSetGenerator {

private:

	std::pair<Eigen::VectorXi, Eigen::VectorXi> bv(SpinInstruction& op, int nQubits) {
		Eigen::VectorXi vx = Eigen::VectorXi::Zero(nQubits);
		Eigen::VectorXi vz = Eigen::VectorXi::Zero(nQubits);

		for (auto term : op.getTerms()) {
			if (term.second == "X") {
				vx(term.first) += 1;
			} else if (term.second == "Z") {
				vz(term.first) += 1;
			} else if (term.second == "Y") {
				vx(term.first) += 1;
				vz(term.first) += 1;
			}
		}

		return std::make_pair(vx, vz);
	};

	int bv_commutator(SpinInstruction& term1, SpinInstruction& term2, int nQubits) {
			auto pair1 = bv(term1, nQubits);
			auto pair2 = bv(term2, nQubits);
			auto scalar = pair1.first.dot(pair2.second) + pair1.second.dot(pair2.first);
			return scalar % 2;
		};

public:

	std::vector<std::vector<int>> getCommutingSet(CompositeSpinInstruction& composite, int n_qubits) {


		std::vector<std::vector<int>> commuting_ops;
		for (int i = 0; i < composite.nInstructions(); i++) {

			auto t_i = std::dynamic_pointer_cast<SpinInstruction>(
					composite.getInstruction(i));
			if (i == 0) {
				commuting_ops.push_back(std::vector<int> { i });
			} else {
				auto comm_ticker = 0;
				for (int j = 0; j < commuting_ops.size(); j++) {
					auto j_op_list = commuting_ops[j];
					int sum = 0;
					for (auto j_op : j_op_list) {
						auto t_jopPtr = std::dynamic_pointer_cast<
								SpinInstruction>(
								composite.getInstruction(j_op));
						sum += bv_commutator(*t_i.get(), *t_jopPtr.get(),
								n_qubits);
					}

					if (sum == 0) {
						commuting_ops[j].push_back(i);
						comm_ticker += 1;
						break;
					}
				}

				if (comm_ticker == 0) {
					commuting_ops.push_back(std::vector<int> { i });
				}
			}
		}

//		for (auto s : commuting_ops) {
//			for (auto x : s) {
//				std::cout << x << " ";
//			}
//			std::cout << "\n";
//		}

		return commuting_ops;
	}

};

}
}

#endif
