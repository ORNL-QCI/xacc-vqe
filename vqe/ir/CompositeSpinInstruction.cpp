#include "SpinInstruction.hpp"

namespace xacc {
namespace vqe {

bool CompositeSpinInstruction::operator==(SpinInstruction &b) {
	if (nInstructions() > 1) {
		return false;
	} else {
		auto casted = std::dynamic_pointer_cast<SpinInstruction>(
				getInstruction(0));
		return casted->operator==(b);
	}
}

bool CompositeSpinInstruction::operator!=(SpinInstruction &b) {
	return !operator==(b);
}

bool CompositeSpinInstruction::operator ==(CompositeSpinInstruction &b) {
	if (nInstructions() != b.nInstructions()) {
		return false;
	}

	for (int i = 0; i < nInstructions(); i++) {
		auto casted1 = std::dynamic_pointer_cast<SpinInstruction>(
				getInstruction(i));
		auto casted2 = std::dynamic_pointer_cast<SpinInstruction>(
				b.getInstruction(i));
		if (casted1 != casted2) {
			return false;
		}
	}

	return true;
}

bool CompositeSpinInstruction::operator!=(CompositeSpinInstruction & b) {
	return !operator==(b);
}

CompositeSpinInstruction CompositeSpinInstruction::operator*(
		SpinInstruction &b) {

	// Multiply a single term by ( sum of terms... )
	CompositeSpinInstruction ret;

	for (auto i : getInstructions()) {
		auto casted = std::dynamic_pointer_cast<SpinInstruction>(i);
		auto newinst = casted->operator*(b);
		ret.addInstruction(std::make_shared<SpinInstruction>(newinst));
	}
	ret.simplify();
	return ret;
}


CompositeSpinInstruction CompositeSpinInstruction::operator*(
		CompositeSpinInstruction &b) {

	CompositeSpinInstruction ret;

	// ( terms... ) * ( terms... )

	for (auto i : getInstructions()) {
		for (auto j : b.getInstructions()) {
			auto casted1 = std::dynamic_pointer_cast<SpinInstruction>(i);
			auto casted2 = std::dynamic_pointer_cast<SpinInstruction>(j);

			auto newinst = casted1->operator*(*casted2.get());

			ret.addInstruction(std::make_shared<SpinInstruction>(newinst));
		}
	}

	ret.simplify();
	return ret;
}

CompositeSpinInstruction CompositeSpinInstruction::operator*(const std::complex<double> &b) {

	// Multiply a single term by ( sum of terms... )
	CompositeSpinInstruction ret;

	for (auto i : getInstructions()) {
		auto casted = std::dynamic_pointer_cast<SpinInstruction>(i);
		auto newinst = casted->operator*(b);
		ret.addInstruction(std::make_shared<SpinInstruction>(newinst));
	}

	ret.simplify();
	return ret;
}

CompositeSpinInstruction CompositeSpinInstruction::operator*(const double &b) {
	return this->operator*(std::complex<double>(b, 0.0));
}

CompositeSpinInstruction& CompositeSpinInstruction::operator*=(const std::complex<double> &b) {
	for (auto i : getInstructions()) {
		auto casted = std::dynamic_pointer_cast<SpinInstruction>(i);
		casted->operator*=(b);
	}
	simplify();
	return *this;
}

CompositeSpinInstruction& CompositeSpinInstruction::operator*=(const double &b) {
	return operator*=(std::complex<double>(b, 0.0));
}

CompositeSpinInstruction CompositeSpinInstruction::operator+(CompositeSpinInstruction& b) {
	CompositeSpinInstruction ret(*this);

	for (auto i : b.getInstructions()) {
		ret.addInstruction(i);
	}

	ret.simplify();

	return ret;
}


void CompositeSpinInstruction::simplify() {
	for (int i = 0; i < nInstructions(); i++) {
		for (int j = 0; j < nInstructions(); j++) {
			if (i < j) {
				auto castedi = std::dynamic_pointer_cast<SpinInstruction>(
						getInstruction(i));
				auto castedj = std::dynamic_pointer_cast<SpinInstruction>(
						getInstruction(j));

				if (castedi->operator==(*castedj.get())) {

					auto size = castedi->bits().size();
					auto coeffi = boost::get<std::complex<double>>(castedi->getParameter(size));
					auto coeffj = boost::get<std::complex<double>>(castedj->getParameter(size));
					auto newParam = InstructionParameter(coeffi + coeffj);

					castedi->setParameter(castedi->bits().size(), newParam);
					auto it = instructions.begin();
					std::advance(it, j);
					instructions.erase(it);
					j--;
				}
			}
		}
	}

	for (int i = 0; i < nInstructions(); i++) {
		auto inst = getInstruction(i);
		auto coeff = boost::get<std::complex<double>>(inst->getParameter(inst->bits().size()));
		if (std::complex<double>(0,0) == coeff) {
			auto it = instructions.begin();
			std::advance(it, i);
			instructions.erase(it);
			i--;
		}
	}
}
}
}


xacc::vqe::CompositeSpinInstruction operator*(double const& scalar,
		xacc::vqe::CompositeSpinInstruction rhs) {
	return rhs * scalar;
}

xacc::vqe::CompositeSpinInstruction operator*(
		std::complex<double> const& scalar,
		xacc::vqe::CompositeSpinInstruction rhs) {
	return rhs * scalar;
}
