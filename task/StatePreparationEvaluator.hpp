/***********************************************************************************
 * Copyright (c) 2017, UT-Battelle
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of the xacc nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Contributors:
 *   Initial API and implementation - Alex McCaskey
 *
 **********************************************************************************/
#ifndef TASK_STATEPREPARATIONEVALUATOR_HPP_
#define TASK_STATEPREPARATIONEVALUATOR_HPP_

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include "exprtk.hpp"
#include <boost/algorithm/string.hpp>
#include "IRProvider.hpp"
#include "XACC.hpp"

namespace xacc {
namespace vqe {

class StatePreparationEvaluator {

protected:

	static constexpr double pi = boost::math::constants::pi<double>();

	using symbol_table_t = exprtk::symbol_table<double>;
	using expression_t = exprtk::expression<double>;
	using parser_t = exprtk::parser<double>;

public:

	static std::shared_ptr<Function> evaluateCircuit(
			const std::shared_ptr<Function> statePrep, const int nParameters,
			const Eigen::VectorXd& x) {
		auto pi = boost::math::constants::pi<double>();
		std::vector<std::string> variableNames;
		for (int i = 0; i < nParameters; i++) {
			variableNames.push_back(
					boost::get<std::string>(statePrep->getParameter(i)));
		}

		auto gateRegistry = xacc::getService<IRProvider>("gate");
		auto evaluatedStatePrep = gateRegistry->createFunction("evaled_"+statePrep->name(), {}, {});

		for (auto inst : statePrep->getInstructions()) {
			if (inst->isParameterized()
					&& inst->getParameter(0).which() == 3) {
				int idx = -1;
				auto expression = boost::get<std::string>(
						inst->getParameter(0));
				for (int i = 0; i < nParameters; i++) {

					if (boost::contains(expression, variableNames[i])) {
						idx = i;
					}
				}

				std::string varName = variableNames[idx];
				double val;
				symbol_table_t symbol_table;
				symbol_table.add_variable(varName, val);
				symbol_table.add_constants();
				expression_t expr;
				expr.register_symbol_table(symbol_table);
				parser_t parser;
				parser.compile(expression, expr);
				val = x(idx);
				auto res = expr.value();
//				if (res < 0.0) {
//					res = 4 * pi + res;
//				}
				InstructionParameter p(res);
				auto updatedInst = gateRegistry->createInstruction(
						inst->name(), inst->bits());
				updatedInst->setParameter(0, p);
				evaluatedStatePrep->addInstruction(updatedInst);
			} else {
				evaluatedStatePrep->addInstruction(inst);
			}
		}

		return evaluatedStatePrep;
	}
};
}
}


#endif /* TASK_STATEPREPARATIONEVALUATOR_HPP_ */
