/*
 * VQEGateFunction.cpp
 *
 *  Created on: Aug 8, 2017
 *      Author: aqw
 */

#include "VQEGateFunction.hpp"
#include <boost/python.hpp>
using namespace boost::python;

namespace xacc {
namespace vqe {
void VQEGateFunction::evaluateVariableParameters(
		std::vector<InstructionParameter> runtimeParameters) {

	Py_Initialize();

	const std::string code = R"code(val = eval(exprStr))code";

	// Retrieve the main module.
	object main = import("__main__");
	object main_namespace = main.attr("__dict__");

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

		boost::trim(oldExprStr);

		main_namespace["exprStr"] = oldExprStr;
		exec(code.c_str(), main_namespace);
		double res = extract<double>(main_namespace["val"]);
		InstructionParameter newParameter(res);

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
						std::string variableExpressionStr = boost::get<
								std::string>(instParam);

						auto runtimeValue = runtimeParameters[i];
						ss1 << runtimeValue;
						boost::replace_all(variableExpressionStr,
								boost::get<std::string>(fParam), ss1.str());

						main_namespace["exprStr"] = variableExpressionStr;
						exec(code.c_str(), main_namespace);
						double res = extract<double>(main_namespace["val"]);
						InstructionParameter oldParameter = inst->getParameter(
								j);
						InstructionParameter newParameter(res);

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

	Py_Finalize();
}
}
}

