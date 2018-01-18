/*
 * StatePreparationEvaluator.hpp
 *
 *  Created on: Nov 7, 2017
 *      Author: aqw
 */

#ifndef TASK_STATEPREPARATIONEVALUATOR_HPP_
#define TASK_STATEPREPARATIONEVALUATOR_HPP_

#include "GateInstruction.hpp"
#include "GateFunction.hpp"
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include "exprtk.hpp"
#include <boost/algorithm/string.hpp>

namespace xacc {
namespace vqe {

using namespace xacc::quantum;

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
		for (int i = 0; i < x.rows(); i++) {
			variableNames.push_back(
					boost::get<std::string>(statePrep->getParameter(i)));
		}

		auto evaluatedStatePrep = std::make_shared<GateFunction>("stateprep");
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
				if (res < 0.0) {
					res = 4 * pi + res;
				}
				InstructionParameter p(res);
				auto updatedInst = GateInstructionRegistry::instance()->create(
						inst->getName(), inst->bits());
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
