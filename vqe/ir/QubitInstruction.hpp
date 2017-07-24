/*
 * QubitInstruction.hpp
 *
 *  Created on: Jul 23, 2017
 *      Author: aqw
 */

#ifndef VQE_IR_QUBITINSTRUCTION_HPP_
#define VQE_IR_QUBITINSTRUCTION_HPP_



#include "Instruction.hpp"

namespace xacc {

namespace vqe {

class QubitInstruction : public Instruction {

public:

	/**
	 * The list of qubits this term operates on
	 */
	std::vector<int> qubits;

	/**
	 * For the ith qubit, the ith InstructionParameter
	 * is a string that is X, Y, or Z. The last
	 * element of this list is the coefficient of this term
	 */
	std::vector<InstructionParameter> parameters;

public:

	/**
	 * The Constructor, takes one qubit
	 * indicating this is a bias value, initialized to 0.0
	 *
	 * @param qbit The bit index
	 */
	QubitInstruction(std::vector<std::pair<int,std::string>> operators) {
		for (auto p : operators) {
			qubits.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(InstructionParameter(1.0));
	}

	QubitInstruction(std::vector<std::pair<int,std::string>> operators, double coeff) {
		for (auto p : operators) {
			qubits.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(InstructionParameter(coeff));
	}

	/**
	 * Return the name of this Instruction
	 *
	 * @return name The name of this Instruction
	 */
	virtual const std::string getName() {
		return "qubit-instruction";
	};

	/**
	 * Persist this Instruction to an assembly-like
	 * string.
	 *
	 * @param bufferVarName The name of the AcceleratorBuffer
	 * @return str The assembly-like string.
	 */
	virtual const std::string toString(const std::string& bufferVarName) {
		std::stringstream ss;
		ss << getParameter(qubits.size()) << " * ";
		for (int i = 0; i < qubits.size(); i++) {
			ss << boost::get<std::string>(getParameter(i)) << qubits[i] << " * ";
		}
		return ss.str().substr(0, ss.str().size() - 2);
	}

	/**
	 * Return the indices of the bits that this Instruction
	 * operates on.
	 *
	 * @return bits The bits this Instruction operates on.
	 */
	virtual const std::vector<int> bits() {
		return qubits;
	}

	/**
	 * Return this Instruction's parameter at the given index.
	 *
	 * @param idx The index of the parameter.
	 * @return param The InstructionParameter at the given index.
	 */
	virtual InstructionParameter getParameter(const int idx) {
		return parameters[idx];
	}

	/**
	 * Return all of this Instruction's parameters.
	 *
	 * @return params This instructions parameters.
	 */
	virtual std::vector<InstructionParameter> getParameters() {
		return parameters;
	}

	/**
	 * Set this Instruction's parameter at the given index.
	 *
	 * @param idx The index of the parameter
	 * @param inst The instruction.
	 */
	virtual void setParameter(const int idx, InstructionParameter& inst) {
		parameters[idx] = inst;
	}

	/**
	 * Return the number of InstructionParameters this Instruction contains.
	 *
	 * @return nInsts The number of instructions.
	 */
	virtual const int nParameters() {
		return parameters.size();
	}

	/**
	 * Return true if this Instruction is parameterized.
	 *
	 * @return parameterized True if this Instruction has parameters.
	 */
	virtual bool isParameterized() {
		return true;
	}
	/**
	 * Returns true if this Instruction is composite,
	 * ie, contains other Instructions.
	 *
	 * @return isComposite True if this is a composite Instruction
	 */
	virtual bool isComposite() { return false; }

	/**
	 * Returns true if this Instruction is enabled
	 *
	 * @return enabled True if this Instruction is enabled.
	 */
	virtual bool isEnabled() { return true; }

	/**
	 * Disable this Instruction
	 */

	virtual void disable() {  }

	/**
	 * Enable this Instruction.
	 */
	virtual void enable() {  }

	void multiplyTerm(QubitInstruction& inst) {

	}

	EMPTY_DEFINE_VISITABLE()
};
}
}


#endif /* VQE_IR_QUBITINSTRUCTION_HPP_ */
