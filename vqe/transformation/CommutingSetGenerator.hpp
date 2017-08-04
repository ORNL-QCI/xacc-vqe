/*
 * CommutingSetGenerator.hpp
 *
 *  Created on: Aug 4, 2017
 *      Author: aqw
 */

#ifndef VQE_TRANSFORMATION_COMMUTINGSETGENERATOR_HPP_
#define VQE_TRANSFORMATION_COMMUTINGSETGENERATOR_HPP_

#include "SpinInstruction.hpp"

namespace xacc {

namespace vqe {
class CommutingSetGenerator {

public:

	std::vector<std::vector<int>> getCommutingSet(CompositeSpinInstruction& instruction) {
		struct IsIIn {
			const int d;
			IsIIn(int n) :
					d(n) {
			}
			bool operator()(std::vector<int> set) const {
				return std::find(set.begin(), set.end(), d) != set.end();
			}
		};
		// Compute commuting sets...
		auto supports = [](CompositeSpinInstruction& i, int idx) -> std::vector<std::pair<int, std::string>> {
			auto spinInst = std::dynamic_pointer_cast<SpinInstruction>(i.getInstruction(idx));
			return spinInst->getTerms();
		};

		auto commutator =
				[](std::vector<std::pair<int, std::string>> support1,
						std::vector<std::pair<int, std::string>> support2) -> bool {
			std::vector<int> iqbits, jqbits;
			std::vector<std::string> ips, jps;

			for (auto t : support1) {
				iqbits.push_back(t.first);
				ips.push_back(t.second);
			}
			for (auto t : support2) {
				jqbits.push_back(t.first);
				jps.push_back(t.second);
			}

			std::vector<int> overlaps;
			for (int i = 0; i < iqbits.size(); i++) {
				auto site = iqbits[i];
				auto itr = std::find(jqbits.begin(), jqbits.end(), site);
				if (itr != jqbits.end()) {
					auto ind2 = std::distance(jqbits.begin(), itr);
					if (ips[i] != jps[ind2]) {
						overlaps.push_back(site);
					}
				}
			}

			return overlaps.size() == 0;
		};

		std::vector<std::vector<int>> commutingSets;
		for (int i = 0; i < instruction.getInstructions().size(); i++) {
			if (i == 0) {
				commutingSets.push_back(std::vector<int>{i});
			}

			for (int j = 0 ; j < i; j++) {
				auto s1 = supports(instruction, i);
				auto s2 = supports(instruction, j);
				if (commutator(s1, s2)) {
					int counter = 0;
					for (auto s : commutingSets) {
						auto itr = std::find(s.begin(), s.end(), j);
						if (itr != s.end()) {
							commutingSets[counter].push_back(i);
							break;
						}
						counter++;
					}
				}
			}



			if (!std::any_of(commutingSets.begin(), commutingSets.end(), IsIIn(i))) {
				commutingSets.push_back(std::vector<int>{i});
			}

		}

		std::cout << "FINAL Commuting term indices:\n";
		for (auto s : commutingSets) {
			std::cout << "{ ";
			for (auto si : s) {
				std::cout << si << " ";
			}
			std::cout << "}";
		}

		return commutingSets;
	}

};

}
}

#endif
