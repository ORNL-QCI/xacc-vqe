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

	FermionInstruction(std::vector<std::pair<int, int>> operators,
			std::string var, std::complex<double> coeff) :
			terms(operators) {
		for (auto p : operators) {
			sites.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(InstructionParameter(coeff));
		parameters.push_back(InstructionParameter(var));
	}

	/**
	 * Return the name of this Instruction
	 *
	 * @return name The name of this Instruction
	 */
	const std::string name() const override {
		return "fermion-instruction";
	};

	const std::string description() const override {
		return "";
	}

    const int nRequiredBits() const override {return 0;}

	void mapBits(std::vector<int> bitMap) override {
	}

	/**
	 * Persist this Instruction to an assembly-like
	 * string.
	 *
	 * @param bufferVarName The name of the AcceleratorBuffer
	 * @return str The assembly-like string.
	 */
	const std::string toString(const std::string& bufferVarName) override {
		std::stringstream ss;
		ss << getParameter(nParameters()-2).toString() << " * ";
		auto variable = getParameter(nParameters()-1).toString();
		if (!variable.empty()) {
			ss << variable << " * ";
		}
		for (int i = 0; i < sites.size(); i++) {
			ss << "a" << sites[i] << ((getParameter(i).as<int>()) ? "\u2020" : "") << " * ";
		}
		return ss.str().substr(0, ss.str().size() - 2);
	}

    const std::string toString() override {
        return toString("");
    }
	/**
	 * Return the indices of the bits that this Instruction
	 * operates on.
	 *
	 * @return bits The bits this Instruction operates on.
	 */
	const std::vector<int> bits() override {
		return sites;
	}

	/**
	 * Return this Instruction's parameter at the given index.
	 *
	 * @param idx The index of the parameter.
	 * @return param The InstructionParameter at the given index.
	 */
	InstructionParameter getParameter(const int idx) const {
		return parameters[idx];
	}

	/**
	 * Return all of this Instruction's parameters.
	 *
	 * @return params This instructions parameters.
	 */
	std::vector<InstructionParameter> getParameters() override {
		return parameters;
	}

	/**
	 * Set this Instruction's parameter at the given index.
	 *
	 * @param idx The index of the parameter
	 * @param inst The instruction.
	 */
	void setParameter(const int idx, InstructionParameter& inst) override {
		parameters[idx] = inst;
	}

	/**
	 * Return the number of InstructionParameters this Instruction contains.
	 *
	 * @return nInsts The number of instructions.
	 */
	const int nParameters() override {
		return parameters.size();
	}

	/**
	 * Return true if this Instruction is parameterized.
	 *
	 * @return parameterized True if this Instruction has parameters.
	 */
	bool isParameterized() override {
		return true;
	}
	/**
	 * Returns true if this Instruction is composite,
	 * ie, contains other Instructions.
	 *
	 * @return isComposite True if this is a composite Instruction
	 */
	bool isComposite() override { return false; }

	/**
	 * Returns true if this Instruction is enabled
	 *
	 * @return enabled True if this Instruction is enabled.
	 */
	bool isEnabled() override { return true; }

	/**
	 * Disable this Instruction
	 */
	void disable() override {  }

	/**
	 * Enable this Instruction.
	 */
	void enable() override {  }

  /**
   * Return true if this Instruction has
   * customizable options.
   *
   * @return hasOptions
   */
  bool hasOptions() override {
      return false;
  }

  /**
   * Set the value of an option with the given name.
   *
   * @param optName The name of the option.
   * @param option The value of the option
   */
  void setOption(const std::string optName,
                 InstructionParameter option) override {
      XACCLogger::instance()->error("setOption not implemented for FermionInst."); 
      return;              
  }
  
  /**
   * Get the value of an option with the given name.
   *
   * @param optName Then name of the option.
   * @return option The value of the option.
   */
  InstructionParameter getOption(const std::string optName) override {
       XACCLogger::instance()->error("getOption not implemented for FermionInst.");  
       return InstructionParameter(0);             
  }

  /**
   * Return all the Instructions options as a map.
   *
   * @return optMap The options map.
   */
  std::map<std::string, InstructionParameter> getOptions() override {
       XACCLogger::instance()->error("getOptions not implemented for FermionInst."); 
       return std::map<std::string,InstructionParameter>();              
  }


	EMPTY_DEFINE_VISITABLE()
};

}

}
#endif
