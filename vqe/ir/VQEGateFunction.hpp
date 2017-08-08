/*
 * VQEGateFunction.hpp
 *
 *  Created on: Aug 8, 2017
 *      Author: aqw
 */

#ifndef VQE_IR_VQEGATEFUNCTION_HPP_
#define VQE_IR_VQEGATEFUNCTION_HPP_

#include "GateFunction.hpp"
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/symbol.h>
#include <symengine/subs.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/parser.h>
#include <symengine/series.h>

#include <boost/algorithm/string.hpp>

using SymEngine::Basic;
using SymEngine::Add;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Symbol;
using SymEngine::symbol;
using SymEngine::umap_basic_num;
using SymEngine::Integer;
using SymEngine::integer;
using SymEngine::multinomial_coefficients;
using SymEngine::RCP;
using SymEngine::rcp_dynamic_cast;
using SymEngine::parse;
using SymEngine::series;
using SymEngine::map_basic_basic;
using SymEngine::subs;

namespace xacc {

namespace vqe {

using namespace xacc::quantum;

/**
 * The VQE Gate Function is a GateFunction that keeps track of the
 * Hamiltonian term leading coefficient, which will be used in
 * the computation of the ground state energy in the wider VQE solve.
 *
 * It also provides an improved evaluateVariableParameters method
 * that takes into account complex variable expressions.
 */
class VQEGateFunction : public GateFunction {

protected:

	std::vector<
			std::pair<std::shared_ptr<Instruction>,
					std::pair<std::pair<int, int>, InstructionParameter>>> cachedInstructions;

public:

	double coefficient = 0.0;

	VQEGateFunction(const std::string& name,
			std::vector<InstructionParameter> varParams, double coeff) :
			GateFunction(name, varParams), coefficient(coeff) {
	}

	VQEGateFunction(GateFunction& gf,
			std::vector<InstructionParameter> varParams, double coeff) :
			 coefficient(coeff), GateFunction(
					gf.getName(), varParams) {
		this->instructions = gf.getInstructions();
	}


	virtual void evaluateVariableParameters(
			std::vector<InstructionParameter> runtimeParameters) {

		std::map<std::string, InstructionParameter> varToValMap;

		assert(runtimeParameters.size() == parameters.size());

		for (auto& foundVarInsts : cachedInstructions) {
			std::stringstream ss, ss1;
			auto inst = foundVarInsts.first;
			auto runtimeParamsIdx = foundVarInsts.second.first.second;
			auto fParamIdx = foundVarInsts.second.first.first;
			auto oldInstParam = foundVarInsts.second.second;

			auto oldExprStr = boost::get<std::string>(oldInstParam);
			auto fParamVarStr = boost::get<std::string>(parameters[fParamIdx]);

			auto runtimeValue = runtimeParameters[fParamIdx];

			ss1 << runtimeValue;
			boost::replace_all(oldExprStr, fParamVarStr, ss1.str());

			auto expr = SymEngine::parse(oldExprStr);
			ss << *expr;
			InstructionParameter newParameter(std::stod(ss.str()));

			inst->setParameter(runtimeParamsIdx, newParameter);

		}

		for (int i = 0; i < parameters.size(); i++) {
			auto fParam = parameters[i];
			for (auto inst : instructions) {
				if (inst->isComposite()) {
					std::dynamic_pointer_cast<Function>(inst)->evaluateVariableParameters(
							runtimeParameters);
				} else if (inst->isParameterized()
						&& inst->getName() != "Measure") {

					for (int j = 0; j < inst->nParameters(); ++j) {
						auto instParam = inst->getParameter(j);

						// If this Instruction's jth parameter is variable (ie is a string)
						// and it contains the current fParam
						if (instParam.which() == 3
								&& boost::contains(
										boost::get<std::string>(instParam),
										boost::get<std::string>(fParam))) {

							std::stringstream ss, ss1;

							// This is a variable
							auto variableExpressionStr = boost::get<std::string>(instParam);

							auto runtimeValue = runtimeParameters[j];
							ss1 << runtimeValue;
							boost::replace_all(variableExpressionStr, boost::get<std::string>(fParam), ss1.str());

							auto expr = SymEngine::parse(variableExpressionStr);
							ss << *expr;

							InstructionParameter oldParameter = inst->getParameter(j);
							InstructionParameter newParameter(std::stod(ss.str()));

							inst->setParameter(j, newParameter);

							cachedInstructions.push_back(
									std::make_pair(inst,
											std::make_pair(
													std::pair<int, int>(i, j),
													oldParameter)));
						}
					}
				}
			}
		}
	}

};

}
}

#endif /* VQE_IR_VQEGATEFUNCTION_HPP_ */
