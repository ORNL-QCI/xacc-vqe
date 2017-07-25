#ifndef VQE_IR_SPININSTRUCTION_HPP_
#define VQE_IR_SPININSTRUCTION_HPP_

#include "Instruction.hpp"
#include "CompositeSpinInstruction.hpp"
#include <map>

namespace xacc {

namespace vqe {

class SpinInstruction: public Instruction {

protected:

	/**
	 * The list of qubits this term operates on
	 */
	std::vector<int> qubits;

	/**
	 * For the ith qubit, the ith InstructionParameter
	 * is a string that is I, X, Y, or Z. The last
	 * element of this list is the coefficient of this term
	 */
	std::vector<InstructionParameter> parameters;

	std::map<std::pair<std::string, std::string>,
			std::pair<std::complex<double>, std::string>> pauliProducts;

public:

	SpinInstruction(const SpinInstruction& i) :
			qubits(i.qubits), parameters(i.parameters), pauliProducts(
					i.pauliProducts) {
	}
	/**
	 * The Constructor, takes one qubit
	 * indicating this is a bias value, initialized to 0.0
	 *
	 * @param qbit The bit index
	 */
	SpinInstruction(std::vector<std::pair<int, std::string>> operators) {
		for (auto p : operators) {
			qubits.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(
				InstructionParameter(std::complex<double>(1.0, 0.0)));
		pauliProducts.insert(
				std::make_pair(std::make_pair("I", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("I", "X"),
						std::make_pair(std::complex<double>(1.0, 0.0), "X")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("X", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "X")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("I", "Y"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Y")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Y", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Y")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Z", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Z")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("I", "Z"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Z")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("X", "X"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Y", "Y"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Z", "Z"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("X", "Y"),
						std::make_pair(std::complex<double>(0.0, 1.0), "Z")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("X", "Z"),
						std::make_pair(std::complex<double>(0.0, -1.0), "Y")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Y", "X"),
						std::make_pair(std::complex<double>(0.0, -1.0), "Z")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Y", "Z"),
						std::make_pair(std::complex<double>(0.0, 1.0), "X")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Z", "X"),
						std::make_pair(std::complex<double>(0.0, 1.0), "Y")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Z", "Y"),
						std::make_pair(std::complex<double>(0.0, -1.0), "X")));
	}

	SpinInstruction(std::vector<std::pair<int, std::string>> operators,
			std::complex<double> coeff) {
		for (auto p : operators) {
			qubits.push_back(p.first);
			parameters.push_back(InstructionParameter(p.second));
		}
		parameters.push_back(InstructionParameter(coeff));
		pauliProducts.insert(
				std::make_pair(std::make_pair("I", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("I", "X"),
						std::make_pair(std::complex<double>(1.0, 0.0), "X")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("X", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "X")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("I", "Y"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Y")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Y", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Y")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Z", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Z")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("I", "Z"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Z")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("X", "X"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Y", "Y"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Z", "Z"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("X", "Y"),
						std::make_pair(std::complex<double>(0.0, 1.0), "Z")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("X", "Z"),
						std::make_pair(std::complex<double>(0.0, -1.0), "Y")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Y", "X"),
						std::make_pair(std::complex<double>(0.0, -1.0), "Z")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Y", "Z"),
						std::make_pair(std::complex<double>(0.0, 1.0), "X")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Z", "X"),
						std::make_pair(std::complex<double>(0.0, 1.0), "Y")));
		pauliProducts.insert(
				std::make_pair(std::make_pair("Z", "Y"),
						std::make_pair(std::complex<double>(0.0, -1.0), "X")));

	}

	/**
	 * Return the name of this Instruction
	 *
	 * @return name The name of this Instruction
	 */
	virtual const std::string getName() {
		return "qubit-instruction";
	}
	;

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
			ss << boost::get<std::string>(getParameter(i)) << qubits[i]
					<< " * ";
		}
		auto r = ss.str().substr(0, ss.str().size() - 2);
		boost::trim(r);
		return r;
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
	virtual bool isComposite() {
		return false;
	}

	/**
	 * Returns true if this Instruction is enabled
	 *
	 * @return enabled True if this Instruction is enabled.
	 */
	virtual bool isEnabled() {
		return true;
	}

	/**
	 * Disable this Instruction
	 */

	virtual void disable() {
	}

	/**
	 * Enable this Instruction.
	 */
	virtual void enable() {
	}

	/**
	 * Return true if the given SpinInstruction is equal to this one
	 *
	 * @param b Other SpinInstruction
	 * @return equal
	 */
	bool operator ==(const SpinInstruction &b) const {
		if (b.qubits.size() != qubits.size()) {
			return false;
		}

		for (int i = 0; i < qubits.size() - 1; i++) {
			if ((qubits[i] != b.qubits[i])
					|| (getParameter(i) != b.getParameter(i))) {
				return false;
			}
		}

		return true;
	}

	/**
	 * Return true if this SpinInstruction is not equal to the provided one
	 * @param b The other SpinInstruction
	 * @return equal
	 */
	bool operator !=(const SpinInstruction &b) const {
		return !operator==(b);
	}

	/**
	 * Multiply this SpinInstruction by the given one, and
	 * return a new SpinInstruction instance.
	 *
	 * @param b The other SpinInstruciton
	 * @return multSpinInstruction
	 */
	SpinInstruction operator*(const SpinInstruction &b) const {

		auto newCoeff = boost::get<std::complex<double>>(
				getParameter(qubits.size()))
				* boost::get<std::complex<double>>(
						b.getParameter(b.qubits.size()));

		std::vector<std::pair<int, std::string>> newTerms;
		for (int i = 0; i < qubits.size(); i++) {
			newTerms.push_back(
					std::make_pair(qubits[i],
							boost::get<std::string>(getParameter(i))));
		}
		for (int i = 0; i < b.qubits.size(); i++) {
			newTerms.push_back(
					std::make_pair(b.qubits[i],
							boost::get<std::string>(b.getParameter(i))));
		}

		return replaceCommonPauliProducts(newTerms, newCoeff);

	}

	/**
	 * Multiply this SpinInstruction by the given scalar complex value
	 * and return the result as a new SpinInstruction
	 *
	 * @param b the value
	 * @return multipliedInst
	 */
	SpinInstruction operator*(const std::complex<double> &b) const {
		SpinInstruction ret(*this);
		InstructionParameter p(
				boost::get<std::complex<double>>(getParameter(qubits.size()))
						* b);
		ret.setParameter(qubits.size(), p);
		return ret;
	}

	/**
	 * Multiply this SpinInstruction by the given scalar double value
	 * and return the result as a new SpinInstruction
	 *
	 * @param b the value
	 * @return multipliedInst
	 */

	SpinInstruction operator*(const double &b) const {
		return operator*(std::complex<double>(b, 0.0));
	}

	/**
	 * Multiply this SpinInstruction by the given scalar double value
	 * and return a reference to this SpinInstruction
	 *
	 * @param b the value
	 * @return multipliedInst
	 */
	SpinInstruction& operator*=(const double &d) {
		return operator*=(std::complex<double>(d, 0.0));
	}

	/**
	 * Multiply this SpinInstruction by the given scalar complex value
	 * and return a reference to this SpinInstruction
	 *
	 * @param b the value
	 * @return multipliedInst
	 */
	SpinInstruction& operator*=(const std::complex<double> &d) {
		InstructionParameter p(
				boost::get<std::complex<double>>(getParameter(qubits.size()))
						* d);
		setParameter(qubits.size(), p);
		return *this;
	}

	/**
	 * Multiply this SpinInstruction by the given CompositeSpinInstruction
	 * and return a new CompositeSpinInstruction.
	 *
	 * @param b the CompositeSpinInstruction
	 * @return multipliedInst
	 */
	CompositeSpinInstruction operator*(CompositeSpinInstruction &b) const {

		// Multiply a single term by ( sum of terms... )
		CompositeSpinInstruction ret;

		for (auto i : b.getInstructions()) {
			auto casted = std::dynamic_pointer_cast<SpinInstruction>(i);
			auto newspinInst = operator*(*casted.get());
			ret.addInstruction(std::make_shared<SpinInstruction>(newspinInst));
		}

		return ret;
	}

	/**
	 * Add the given SpinInstruction to this SpinInstruction,
	 * creating a new CompositeSpinInstruction.
	 *
	 * @param b The instruction to add to this one
	 * @return addCompInst The sum result.
	 */
	CompositeSpinInstruction operator+(const SpinInstruction &b) const {
		CompositeSpinInstruction ret;

		if (operator==(b)) {
			auto newCoeff = boost::get<std::complex<double>>(
					getParameter(qubits.size()))
					+ boost::get<std::complex<double>>(
							b.getParameter(qubits.size()));
			auto ptr = std::make_shared<SpinInstruction>(*this);
			InstructionParameter newparam(newCoeff);
			ptr->setParameter(qubits.size(), newparam);
			ret.addInstruction(ptr);
		} else {
			ret.addInstruction(std::make_shared<SpinInstruction>(*this));
			ret.addInstruction(std::make_shared<SpinInstruction>(b));
		}
		return ret;
	}

	EMPTY_DEFINE_VISITABLE()

private:

	SpinInstruction replaceCommonPauliProducts(SpinInstruction& si) {
		auto newCoeff = boost::get<std::complex<double>>(
				getParameter(qubits.size()));
		std::vector<std::pair<int, std::string>> newTerms;
		for (int i = 0; i < si.qubits.size(); i++) {
			newTerms.push_back(
					std::make_pair(qubits[i],
							boost::get<std::string>(getParameter(i))));
		}

		return replaceCommonPauliProducts(newTerms, newCoeff);
	}

	SpinInstruction replaceCommonPauliProducts(
			std::vector<std::pair<int, std::string>>& newTerms,
			std::complex<double> newCoeff) const {

		for (int i = 0; i < newTerms.size() - 1; i++) {
			auto qubit1 = newTerms[i].first;
			auto qubit2 = newTerms[i + 1].first;
			if (qubit1 == qubit2) {
				auto gate1 = newTerms[i].second;
				auto gate2 = newTerms[i + 1].second;
				auto gPair = std::make_pair(gate1, gate2);
				if (pauliProducts.find(gPair) != pauliProducts.end()) {
					newTerms.erase(newTerms.begin() + i + 1);
					auto coeff = pauliProducts.at(gPair).first;
					newTerms.at(i) = std::make_pair(qubit1,
							pauliProducts.at(gPair).second);
					newCoeff *= coeff;
				}
			}
		}

		return SpinInstruction(newTerms, newCoeff);
	}
};
}
}

/**
 * Operator to enable lhs mutliplcation by scalar
 * @param scalar
 * @param rhs
 * @return
 */
xacc::vqe::SpinInstruction operator*(double const& scalar,
		xacc::vqe::SpinInstruction rhs) {
	return rhs *= scalar;
}

xacc::vqe::SpinInstruction operator*(std::complex<double> const& scalar,
		xacc::vqe::SpinInstruction rhs) {
	return rhs *= scalar;
}

#endif /* VQE_IR_SPININSTRUCTION_HPP_ */
