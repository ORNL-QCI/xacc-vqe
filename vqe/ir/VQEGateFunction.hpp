/*
 * VQEGateFunction.hpp
 *
 *  Created on: Aug 8, 2017
 *      Author: aqw
 */

#ifndef VQE_IR_VQEGATEFUNCTION_HPP_
#define VQE_IR_VQEGATEFUNCTION_HPP_

#include "GateFunction.hpp"
#include <boost/algorithm/string.hpp>

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
class VQEGateFunction: public GateFunction {

protected:

	std::vector<
			std::pair<std::shared_ptr<Instruction>,
					std::pair<std::pair<int, int>, InstructionParameter>>> cachedInstructions;

public:

	double coefficient = 0.0;

	bool isIdentityOperator = false;

	VQEGateFunction(const std::string& name,
			std::vector<InstructionParameter> varParams, double coeff) :
			GateFunction(name, varParams), coefficient(coeff) {
	}

	virtual void evaluateVariableParameters(
			std::vector<InstructionParameter> runtimeParameters);

	virtual std::shared_ptr<Function> evaluate(
			std::vector<InstructionParameter> parameters) {
	}
};

}
}

#endif /* VQE_IR_VQEGATEFUNCTION_HPP_ */
