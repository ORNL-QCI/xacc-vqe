/*
 * CompositeSpinInstruction.hpp
 *
 *  Created on: Jul 25, 2017
 *      Author: aqw
 */

#ifndef VQE_IR_COMPOSITESPININSTRUCTION_HPP_
#define VQE_IR_COMPOSITESPININSTRUCTION_HPP_


#include "Function.hpp"
#include <boost/algorithm/string.hpp>

namespace xacc {

namespace vqe {

class CompositeSpinInstruction : public Function {

protected:

	std::list<InstPtr> instructions;

public:

	/**
	 * The constructor, takes the function unique id and its name.
	 *
	 * @param id
	 * @param name
	 */
	CompositeSpinInstruction() {
	}

	virtual const int nInstructions() {
		return instructions.size();
	}

	virtual InstPtr getInstruction(const int idx) {
		if (instructions.size() > idx) {
			return *std::next(instructions.begin(), idx);
		} else {
			XACCError("Invalid instruction index.");
		}
	}

	virtual std::list<InstPtr> getInstructions() {
		return instructions;
	}

	virtual void removeInstruction(const int idx) {
		instructions.remove(getInstruction(idx));
	}

	/**
	 * Add an instruction to this quantum
	 * intermediate representation.
	 *
	 * @param instruction
	 */
	virtual void addInstruction(InstPtr instruction) {
		instructions.push_back(instruction);
	}

	/**
	 * Replace the given current quantum instruction
	 * with the new replacingInst quantum Instruction.
	 *
	 * @param currentInst
	 * @param replacingInst
	 */
	virtual void replaceInstruction(const int idx, InstPtr replacingInst) {
		std::replace(instructions.begin(), instructions.end(),
				getInstruction(idx), replacingInst);
	}

	virtual void insertInstruction(const int idx, InstPtr newInst) {
		auto iter = std::next(instructions.begin(), idx);
		instructions.insert(iter, newInst);
	}

	/**
	 * Return the name of this function
	 * @return
	 */
	virtual const std::string getName() {
		return "";
	}

	/**
	 * Return the qubits this function acts on.
	 * @return
	 */
	virtual const std::vector<int> bits() {
		return std::vector<int> { };
	}

	/**
	 * Return an assembly-like string representation for this function .
	 * @param bufferVarName
	 * @return
	 */
	virtual const std::string toString(const std::string& bufferVarName) {
		std::stringstream ss;
		for (auto i : instructions) {
			ss << i->toString("") << " + ";
		}
		auto r = ss.str().substr(0, ss.str().length() - 2);
		boost::trim(r);
		return r;
	}

	virtual InstructionParameter getParameter(const int idx) const {
		XACCError("CompositeSpinInstruction does not contain runtime parameters.");
	}

	virtual void setParameter(const int idx, InstructionParameter& p) {
		XACCError("CompositeSpinInstruction does not contain runtime parameters.");
	}

	virtual std::vector<InstructionParameter> getParameters() {
		XACCError("CompositeSpinInstruction does not contain runtime parameters.");
	}

	virtual bool isParameterized() {
		return false;
	}

	virtual const int nParameters() {
		return 0;
	}

	virtual void evaluateVariableParameters(std::vector<InstructionParameter> parameters) {
		XACCError("CompositeSpinInstruction does not contain runtime parameters.");
	}


	EMPTY_DEFINE_VISITABLE()

};

}

}
#endif /* VQE_IR_COMPOSITESPININSTRUCTION_HPP_ */
