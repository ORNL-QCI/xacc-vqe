#include "SpinInstruction.hpp"

xacc::vqe::SpinInstruction operator*(double const& scalar,
		xacc::vqe::SpinInstruction rhs) {
	return rhs *= scalar;
}

xacc::vqe::SpinInstruction operator*(std::complex<double> const& scalar,
		xacc::vqe::SpinInstruction rhs) {
	return rhs *= scalar;
}



