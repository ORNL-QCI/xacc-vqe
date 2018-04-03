/***********************************************************************************
 * Copyright (c) 2018, UT-Battelle
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
#ifndef VQE_IR_SPININSTRUCTION_HPP_
#define VQE_IR_SPININSTRUCTION_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <unsupported/Eigen/KroneckerProduct>
#include "XACC.hpp"

// Putting this here due to clang error
// not able to find operator!= from operators.hpp
namespace xacc {
namespace vqe {
class PauliOperator;
}
}
bool operator==(const xacc::vqe::PauliOperator& lhs,
	const xacc::vqe::PauliOperator& rhs);

#include "operators.hpp"

namespace xacc {

namespace vqe {

// A Term can be a coefficient, a variable coefficient, and the terms themselves
using TermTuple = std::tuple<std::complex<double>, std::string, std::map<int, std::string>>;
using c = std::complex<double>;
using ActionResult = std::pair<std::string, c>;
enum ActionType {Bra, Ket};
class Triplet : std::tuple<std::uint64_t, std::uint64_t, std::complex<double>> {
public:
	Triplet(std::uint64_t r, std::uint64_t c, std::complex<double> coeff) {
		std::get<0>(*this) = r;
		std::get<1>(*this) = c;
		std::get<2>(*this) = coeff;
	}
	const std::uint64_t row() {return std::get<0>(*this);}
	const std::uint64_t col() {return std::get<1>(*this);}
	const std::complex<double> coeff() {return std::get<2>(*this);}
};

class Term: public TermTuple,
		public tao::operators::commutative_multipliable<Term>,
		public tao::operators::equality_comparable<Term> {

protected:

	static std::map<std::string, std::pair<c, std::string>> create_map() {
		static std::map<std::string, std::pair<c, std::string>> m;
		m.insert( { "II", { c(1.0, 0.0), "I" } });
		m.insert( { "IX", { c(1.0, 0.0), "X" } });
		m.insert( { "XI", { c(1.0, 0.0), "X" } });
		m.insert( { "IY", { c(1.0, 0.0), "Y" } });
		m.insert( { "YI", { c(1.0, 0.0), "Y" } });
		m.insert( { "ZI", { c(1.0, 0.0), "Z" } });
		m.insert( { "IZ", { c(1.0, 0.0), "Z" } });
		m.insert( { "XX", { c(1.0, 0.0), "I" } });
		m.insert( { "YY", { c(1.0, 0.0), "I" } });
		m.insert( { "ZZ", { c(1.0, 0.0), "I" } });
		m.insert( { "XY", { c(0.0, 1.0), "Z" } });
		m.insert( { "XZ", { c(0.0, -1.0), "Y" } });
		m.insert( { "YX", { c(0.0, -1.0), "Z" } });
		m.insert( { "YZ", { c(0.0, 1.0), "X" } });
		m.insert( { "ZX", { c(0.0, 1.0), "Y" } });
		m.insert( { "ZY", { c(0.0, -1.0), "X" } });
		return m;
	}

	static const std::map<std::string, std::pair<c, std::string>> pauliProducts;

public:

	Term() {
		std::get<0>(*this) = std::complex<double>(0, 0);
		std::get<1>(*this) = "";
		std::get<2>(*this) = { };
	}

	Term(const Term& t) {
		std::get<0>(*this) = std::get<0>(t);
		std::get<1>(*this) = std::get<1>(t);
		std::get<2>(*this) = std::get<2>(t);
	}

	Term(std::complex<double> c) {
		std::get<0>(*this) = c;
		std::get<1>(*this) = "";
		std::get<2>(*this) = { };
	}

	Term(double c) {
		std::get<0>(*this) = std::complex<double>(c, 0);
		std::get<1>(*this) = "";
		std::get<2>(*this) = { };
	}

	Term(std::complex<double> c, std::map<int, std::string> ops) {
		std::get<0>(*this) = c;
		std::get<1>(*this) = "";
		std::get<2>(*this) = ops;
	}

	Term(std::string var) {
		std::get<0>(*this) = std::complex<double>(1,0);
		std::get<1>(*this) = var;
		std::get<2>(*this) = {};
	}

	Term(std::complex<double> c, std::string var) {
		std::get<0>(*this) = c;
		std::get<1>(*this) = var;
		std::get<2>(*this) = { };
	}

	Term(std::string var, std::map<int, std::string> ops) {
		std::get<0>(*this) = std::complex<double>(1,0);
		std::get<1>(*this) = var;
		std::get<2>(*this) = ops;
	}

	Term(std::complex<double> c, std::string var, std::map<int, std::string> ops) {
		std::get<0>(*this) = c;
		std::get<1>(*this) = var;
		std::get<2>(*this) = ops;
	}

	Term(std::map<int, std::string> ops) {
		std::get<0>(*this) = std::complex<double>(1,0);
		std::get<1>(*this) = "";
		std::get<2>(*this) = ops;
	}

	static const std::string id(const std::map<int, std::string>& ops, const std::string& var = "") {
		std::string s;
		s = var;
		for (auto& t : ops) {
			if (t.second != "I") {
				s += t.second + std::to_string(t.first);
			}
		}

		if (s.empty()) {
			return "I";
		}

		return s;
	}

	const std::string id() const {
		std::string s;
		s = std::get<1>(*this);
		for (auto& t : std::get<2>(*this)) {
			if (t.second != "I") {
				s += t.second + std::to_string(t.first);
			}
		}

		if (s.empty()) {
			return "I";
		}

		return s;
	}

	std::map<int, std::string>& ops() {
		return std::get<2>(*this);
	}

	bool isIdentity() {
		if (ops().empty()) {
			return true;
		} else {
			return false;
		}
	}

	std::complex<double>& coeff() {
		return std::get<0>(*this);
	}

	std::string& var() {
		return std::get<1>(*this);
	}

	Term& operator*=( const Term& v ) noexcept;

	bool operator==( const Term& v ) noexcept {
		return std::get<1>(*this) == std::get<1>(v) && ops() == std::get<2>(v);
	}

	std::vector<Triplet> getSparseMatrixElements(const int nQubits);

	ActionResult action(const std::string& bitString, ActionType type);

};

class PauliOperator: public tao::operators::commutative_ring<PauliOperator>,
		public tao::operators::equality_comparable<PauliOperator>,
		public tao::operators::commutative_multipliable<PauliOperator, double>,
		public tao::operators::commutative_multipliable<PauliOperator,
				std::complex<double>> {
protected:

	std::unordered_map<std::string, Term> terms;

public:

	std::unordered_map<std::string, Term>::iterator begin() {
		return terms.begin();
	}
	std::unordered_map<std::string, Term>::iterator end() {
		return terms.end();
	}

	PauliOperator();
	PauliOperator(std::complex<double> c);
	PauliOperator(double c);
	PauliOperator(std::string var);
	PauliOperator(std::complex<double> c, std::string var);
	PauliOperator(const PauliOperator& i);
	PauliOperator(std::map<int, std::string> operators);
	PauliOperator(std::map<int, std::string> operators, std::string var);
	PauliOperator(std::map<int, std::string> operators,
			std::complex<double> coeff);
	PauliOperator(std::map<int, std::string> operators,
			double coeff);
	PauliOperator(std::map<int, std::string> operators,
			std::complex<double> coeff, std::string var);

	const std::vector<std::pair<std::string, std::complex<double>>> computeActionOnKet(
			const std::string& bitString);
	const std::vector<std::pair<std::string, std::complex<double>>> computeActionOnBra(
			const std::string& bitString);

	const int nTerms();

	const std::string toString();

	void clear();

	std::unordered_map<std::string, Term> getTerms() const {
		return terms;
	}

	std::vector<Triplet> getSparseMatrixElements();
	std::shared_ptr<IR> toXACCIR();
	void fromXACCIR(std::shared_ptr<IR> ir);
	PauliOperator eval(const std::map<std::string, std::complex<double>> varToValMap);
	bool isClose(PauliOperator& other);

	PauliOperator& operator+=( const PauliOperator& v ) noexcept;
	PauliOperator& operator-=( const PauliOperator& v ) noexcept;
	PauliOperator& operator*=( const PauliOperator& v ) noexcept;
	bool operator==( const PauliOperator& v ) noexcept;
	PauliOperator& operator*=( const double v ) noexcept;
	PauliOperator& operator*=( const std::complex<double> v ) noexcept;
};
}
}

#endif
