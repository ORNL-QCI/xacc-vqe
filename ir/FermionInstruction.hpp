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

/**
 * The FermionInstruction is a realization of the Instruction
 * interface that models a fermionic creation/annihilation term
 * in a second-quantized Hamiltonian. A FermionInstruction keeps
 * track of a vector of operators, with each operator describing
 * the site it operates on and whether it is a creation or
 * annihilation operator on that site. FermionInstructions
 * also keep track of a coefficient multiplying the term.
 *
 */
class FermionInstruction : public Instruction {

public:

	/**
	 * The individual operator sites. For example,
	 * a^3 a0 - > (3, 0)
	 */
	std::vector<int> sites;
	std::vector<std::pair<int, int>> terms;

	/**
	 * A vector where the ith element represents
	 * a 0 or 1 (annihilation or creation) for the ith
	 * site in the sites vector. The sites.size() element
	 * of this vector represents the term coefficient.
	 */
	std::vector<InstructionParameter> parameters;

public:

	/**
	 * The constructor, takes the (site, creation/annihilation) pairs
	 * that describe this term. Initializes coefficient to 1.0
	 *
	 * @param operators Pairs describing site and if creation or annihilation
	 */
	FermionInstruction(std::vector<std::pair<int, int>> operators) :
			terms(operators) {//, coefficient(std::complex<double>(1.0)) {
		for (auto p : operators) {
			sites.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(InstructionParameter(std::complex<double>(1,0)));
		parameters.push_back(InstructionParameter(std::string("")));

	}

	/**
	 * The constructor, takes the (site, creation/annihilation) pairs
	 * that describe this term and term coefficient
	 *
	 * @param operators Pairs describing site and if creation or annihilation
	 * @param coeff The term coefficient
	 */
	FermionInstruction(std::vector<std::pair<int, int>> operators,
			std::complex<double> coeff) :
			 terms(operators) {
		for (auto p : operators) {
			sites.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(InstructionParameter(coeff));
		parameters.push_back(InstructionParameter(std::string("")));
	}

	FermionInstruction(std::vector<std::pair<int, int>> operators,
			std::string var) :
			terms(operators) {
		for (auto p : operators) {
			sites.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(InstructionParameter(std::complex<double>(1,0)));
		parameters.push_back(InstructionParameter(var));
	}

	/**
	 * Return the name of this Instruction
	 *
	 * @return name The name of this Instruction
	 */
	virtual const std::string name() const {
		return "fermion-instruction";
	};

	virtual const std::string description() const {
		return "";
	}

	virtual void mapBits(std::vector<int> bitMap) {
	}

	/**
	 * Persist this Instruction to an assembly-like
	 * string.
	 *
	 * @param bufferVarName The name of the AcceleratorBuffer
	 * @return str The assembly-like string.
	 */
	virtual const std::string toString(const std::string& bufferVarName) {
		std::stringstream ss;
		ss << getParameter(nParameters()-2) << " * ";
		auto variable = getParameter(nParameters()-1);
		if (!variable.empty()) {
			ss << variable << " * ";
		}
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
	virtual InstructionParameter getParameter(const int idx) const {
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

	virtual const std::string getTag() {
		return "";
	}

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
