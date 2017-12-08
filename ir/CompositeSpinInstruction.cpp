#include "SpinInstruction.hpp"
#include <boost/mpi.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

namespace xacc {
namespace vqe {

using PersistedTriplet = std::pair<std::pair<int,int>, std::complex<double>>;

using Triplet = Eigen::Triplet<std::complex<double>>;
using SparseMat = Eigen::SparseMatrix<std::complex<double>>;

struct add_triplets {
	std::vector<PersistedTriplet> operator()( std::vector<PersistedTriplet> a, std::vector<PersistedTriplet> b) {
		std::move(b.begin(), b.end(), std::back_inserter(a));
		return a;
	}
};

CompositeSpinInstruction::CompositeSpinInstruction() {
}

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
//	ret.simplify();
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

//	ret.simplify();
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

//	ret.simplify();
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
//	simplify();
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

//	ret.simplify();

	return ret;
}

CompositeSpinInstruction CompositeSpinInstruction::operator+(SpinInstruction& b) {
	CompositeSpinInstruction ret(*this);
	ret.addInstruction(std::make_shared<SpinInstruction>(b));
//	ret.simplify();
	return ret;
}

void CompositeSpinInstruction::simplify() {

	std::unordered_map<std::string, std::shared_ptr<SpinInstruction>> map;
	for (int i = 0; i < nInstructions(); ++i) {
		auto inst = std::dynamic_pointer_cast<SpinInstruction>(
								getInstruction(i));
		auto terms = inst->getTerms();
		std::string key;
		for (auto t : terms) {
			if (t.second != "I") {
				key += t.second + std::to_string(t.first);
			} else {
				key = "I";
			}
		}

		auto search = map.find(key);
		if (search != map.end()) {
			search->second->coefficient += inst->coefficient;
		} else {
			map[key] = inst;
		}
	}

	instructions.clear();

	for (auto& kv : map) {
		if (std::fabs(std::real(kv.second->coefficient)) < 1e-12) kv.second->coefficient.real(0.0);
		if (std::fabs(std::imag(kv.second->coefficient)) < 1e-12) kv.second->coefficient.imag(0.0);

		if (kv.second->coefficient != std::complex<double>(0.0,0.0)) {
			addInstruction(kv.second);
		}
	}
}

Eigen::MatrixXcd CompositeSpinInstruction::toMatrix(const int nQubits) {
	int dim = std::pow(2,nQubits);
	Eigen::MatrixXcd hamiltonian = Eigen::MatrixXcd::Zero(dim, dim);
	for (auto instruction : getInstructions()) {
		hamiltonian +=
				std::dynamic_pointer_cast<SpinInstruction>(instruction)->toMatrix(
						nQubits);
	}
	return hamiltonian;
}

Eigen::SparseMatrix<std::complex<double>> CompositeSpinInstruction::toSparseMatrix(const int nQubits) {
	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	Eigen::SparseMatrix<std::complex<double>> ham (dim, dim);
	boost::mpi::communicator world;
	for (int i = 0; i < nInstructions(); i++) {
		auto inst = getInstruction(i);
		ham += std::dynamic_pointer_cast<SpinInstruction>(inst)->toSparseMatrix(
				nQubits);
	}

	return ham;
}

Eigen::SparseMatrix<double> CompositeSpinInstruction::toSparseRealMatrix(const int nQubits) {
	Eigen::SparseMatrix<double> ham = std::dynamic_pointer_cast<
			SpinInstruction>(getInstruction(0))->toSparseRealMatrix(nQubits);
	for (int i = 1; i < nInstructions(); i++) {
		auto inst = getInstruction(i);
		ham += std::dynamic_pointer_cast<SpinInstruction>(inst)->toSparseRealMatrix(
				nQubits);
	}
	return ham;
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
