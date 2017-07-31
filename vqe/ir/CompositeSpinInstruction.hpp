/***********************************************************************************
 * Copyright (c) 2016, UT-Battelle
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
#ifndef VQE_IR_COMPOSITESPININSTRUCTION_HPP_
#define VQE_IR_COMPOSITESPININSTRUCTION_HPP_

#include "Function.hpp"
#include <boost/algorithm/string.hpp>

namespace xacc {

namespace vqe {

class SpinInstruction;

/**
 * The CompositeSpinInstruction is an XACC Function realization
 * that encapsulates a sum of SpinInstructions.
 */
class CompositeSpinInstruction: public Function {

protected:

	/**
	 * List of SpinInstructions this CompositeSpinInstruction
	 * is composed of.
	 */
	std::list<InstPtr> instructions;

public:

	/**
	 * The nullary constructor.
	 */
	CompositeSpinInstruction() {}

	/**
	 * The Copy constructor
	 * @param i
	 */
	CompositeSpinInstruction(const CompositeSpinInstruction& i) :
			instructions(i.instructions) {
	}

	/**
	 * Return the number of SpinInstructions in this sum.
	 *
	 * @return
	 */
	virtual const int nInstructions() {
		return instructions.size();
	}

	/**
	 * Return the SpinInstruction at the given index.
	 *
	 * @param idx The index of the SpinInstruction
	 * @return instruction
	 */
	virtual InstPtr getInstruction(const int idx) {
		if (instructions.size() > idx) {
			return *std::next(instructions.begin(), idx);
		} else {
			XACCError("Invalid instruction index.");
		}
	}

	/**
	 * Return all SpinInstructions
	 *
	 * @return instructions
	 */
	virtual std::list<InstPtr> getInstructions() {
		return instructions;
	}

	/**
	 * Remove the instruction at the given index.
	 *
	 * @param idx The index of the instruction to remove.
	 */
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

	/**
	 * Insert the given SpinInstruction at the given index.
	 *
	 * @param idx The index of the new instruction
	 * @param newInst
	 */
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

		std::string r = "";
		if (instructions.empty()) {
			r = "(0,0)";
		} else {
			for (auto i : instructions) {
				ss << i->toString("") << " + ";
			}
			r = ss.str().substr(0, ss.str().length() - 2);
		}
		boost::trim(r);
		return r;
	}

	/**
	 * CompositeSpinInstructions do not contain parameters. This
	 * method is not implmeneted. An error will be thrown if invoked.
	 *
	 * @param idx
	 * @return
	 */
	virtual InstructionParameter getParameter(const int idx) const {
		XACCError(
				"CompositeSpinInstruction does not contain runtime parameters.");
	}

	/**
	 * CompositeSpinInstructions do not contain parameters. This
	 * method is not implmeneted. An error will be thrown if invoked.
	 *
	 * @param idx
	 * @return
	 */
	virtual void setParameter(const int idx, InstructionParameter& p) {
		XACCError(
				"CompositeSpinInstruction does not contain runtime parameters.");
	}

	/**
	 * CompositeSpinInstructions do not contain parameters. This
	 * method is not implmeneted. An error will be thrown if invoked.
	 *
	 * @param idx
	 * @return
	 */
	virtual std::vector<InstructionParameter> getParameters() {
		XACCError(
				"CompositeSpinInstruction does not contain runtime parameters.");
	}

	/**
	 * CompositeSpinInstruction is not parameterized.
	 *
	 * @return
	 */
	virtual bool isParameterized() {
		return false;
	}

	/**
	 * There are 0 parameters.
	 * @return
	 */
	virtual const int nParameters() {
		return 0;
	}

	/**
	 * CompositeSpinInstructions do not contain parameters. This
	 * method is not implmeneted. An error will be thrown if invoked.
	 *
	 * @param idx
	 * @return
	 */
	virtual void evaluateVariableParameters(
			std::vector<InstructionParameter> parameters) {
		XACCError(
				"CompositeSpinInstruction does not contain runtime parameters.");
	}

	/**
	 * Overloaded operator to check equality between this CompositeSpinInstruction
	 * and SpinInstructions
	 * @param b the spin instruction
	 * @return equal
	 */
	bool operator==(SpinInstruction &b);

	/**
	 * Overloaded operator to check inequality between this CompositeSpinInstruction
	 * and SpinInstructions
	 * @param b the spin instruction
	 * @return equal
	 */
	bool operator!=(SpinInstruction &b);

	/**
	 * Overloaded operator to check equality between this CompositeSpinInstruction
	 * and CompositeSpinInstructions
	 * @param b the composite spin instruction
	 * @return equal
	 */
	bool operator ==(CompositeSpinInstruction &b);

	/**
	 * Overloaded operator to check inequality between this CompositeSpinInstruction
	 * and CompositeSpinInstructions
	 * @param b the composite spin instruction
	 * @return equal
	 */
	bool operator!=(CompositeSpinInstruction & b);

	/**
	 * Multiply this CompositeSpinInstruction by the given SpinInstruction, and
	 * return a new CompositeSpinInstruction instance.
	 *
	 * @param b The other SpinInstruciton
	 * @return multSpinInstruction
	 */
	CompositeSpinInstruction operator*(SpinInstruction &b);

	/**
	 * Multiply this CompositeSpinInstruction by the given CompositeSpinInstruction, and
	 * return a new CompositeSpinInstruction instance.
	 *
	 * @param b The other CompositeSpinInstruciton
	 * @return multSpinInstruction
	 */
	CompositeSpinInstruction operator*(CompositeSpinInstruction &b);

	/**
	 * Multiply this CompositeSpinInstruction by the given scalar complex value
	 * and return the result as a new CompositeSpinInstruction
	 *
	 * @param b the value
	 * @return multipliedInst
	 */
	CompositeSpinInstruction operator*(const std::complex<double> &b);

	/**
	 * Multiply this CompositeSpinInstruction by the given scalar double value
	 * and return the result as a new CompositeSpinInstruction
	 *
	 * @param b the value
	 * @return multipliedInst
	 */
	CompositeSpinInstruction operator*(const double &b);

	/**
	 * Multiply this CompositeSpinInstruction by the given scalar complex value
	 * and return a reference to this CompositeSpinInstruction
	 *
	 * @param b the value
	 * @return multipliedInst
	 */
	CompositeSpinInstruction& operator*=(const std::complex<double> &b);

	/**
	 * Multiply this CompositeSpinInstruction by the given scalar double value
	 * and return a reference to this CompositeSpinInstruction
	 *
	 * @param b the value
	 * @return multipliedInst
	 */
	CompositeSpinInstruction& operator*=(const double &b);

	/**
	 * Add the given CompositeSpinInstruction to this CompositeSpinInstruction,
	 * creating a new CompositeSpinInstruction.
	 *
	 * @param b The instruction to add to this one
	 * @return addCompInst The sum result.
	 */
	CompositeSpinInstruction operator+(CompositeSpinInstruction& b);

	/**
	 * Add the given SpinInstruction to this CompositeSpinInstruction,
	 * creating a new CompositeSpinInstruction.
	 *
	 * @param b The instruction to add to this one
	 * @return addCompInst The sum result.
	 */
	CompositeSpinInstruction operator+(SpinInstruction& b);

	EMPTY_DEFINE_VISITABLE()

private:

	/**
	 * This method simplifies the CompositeSpinInstruction by
	 * adding like terms.
	 */
	void simplify();
};

}

}

xacc::vqe::CompositeSpinInstruction operator*(double const& scalar,
		xacc::vqe::CompositeSpinInstruction rhs);
xacc::vqe::CompositeSpinInstruction operator*(
		std::complex<double> const& scalar,
		xacc::vqe::CompositeSpinInstruction rhs);


#endif
