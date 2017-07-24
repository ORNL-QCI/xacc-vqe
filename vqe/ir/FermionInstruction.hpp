/*
 * FermionInstruction.hpp
 *
 *  Created on: Jul 22, 2017
 *      Author: aqw
 */

#ifndef VQE_IR_FERMIONINSTRUCTION_HPP_
#define VQE_IR_FERMIONINSTRUCTION_HPP_


#include "Instruction.hpp"

namespace xacc {

namespace vqe {

class FermionInstruction : public Instruction {

public:

	/**
	 */
	std::vector<int> sites;

	/**
	 */
	std::vector<InstructionParameter> parameters;

public:

	/**
	 * The Constructor, takes one qubit
	 * indicating this is a bias value, initialized to 0.0
	 *
	 * @param qbit The bit index
	 */
	FermionInstruction(std::vector<std::pair<int,int>> operators) {
		for (auto p : operators) {
			sites.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(InstructionParameter(1.0));
	}

	FermionInstruction(std::vector<std::pair<int,int>> operators, double coeff) {
		for (auto p : operators) {
			sites.push_back(p.first);
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
		return "fermion-instruction";
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
		ss << getParameter(sites.size()) << " * ";
		for (int i = 0; i < sites.size(); i++) {
			ss << "a" << sites[i] << (boost::get<int>(getParameter(i)) ? "\u2020" : "") << " * ";
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
		return sites;
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

	EMPTY_DEFINE_VISITABLE()
};

}

}
#endif
