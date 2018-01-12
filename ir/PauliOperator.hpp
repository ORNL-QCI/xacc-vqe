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
#ifndef VQE_IR_SPININSTRUCTION_HPP_
#define VQE_IR_SPININSTRUCTION_HPP_

#include "CommonPauliProducts.hpp"
#include <map>
#include <unsupported/Eigen/KroneckerProduct>
#include "XACC.hpp"
#include "operators.hpp"

namespace xacc {

namespace vqe {

// A Term can be a coefficient, a variable coefficient, and the terms themselves
using TermTuple = std::tuple<std::complex<double>, std::string, std::map<int, std::string>>;
using PauliProductsType = std::unordered_map<std::string, std::pair<std::complex<double>, std::string>>;
using c = std::complex<double>;

class Term: public TermTuple,
		public tao::operators::commutative_multipliable<Term> {

protected:

	CommonPauliProducts pauliProducts;

public:

	Term() {
		std::get<0>(*this) = std::complex<double>(0, 0);
		std::get<1>(*this) = "";
		std::get<2>(*this) = { };
	}

	Term(std::complex<double> c) {
		std::get<0>(*this) = c;
		std::get<1>(*this) = "";
		std::get<2>(*this) = { {0,"I"}};
	}

	Term(double c) {
		std::get<0>(*this) = std::complex<double>(c, 0);
		std::get<1>(*this) = "";
		std::get<2>(*this) = { {0,"I"}};
	}

	Term(std::complex<double> c, std::map<int, std::string> ops) {
		std::get<0>(*this) = c;
		std::get<1>(*this) = "";
		std::get<2>(*this) = ops;
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

	const std::string id() const {
		std::stringstream s;
		s << std::get<1>(*this);
		for (auto& t : std::get<2>(*this)) {
			if (t.second != "I") {
				s << t.second << t.first;
			}
		}

		if (s.str().empty()) {
			return "I";
		}

		return s.str();
	}

	std::map<int, std::string>& ops() {
		return std::get<2>(*this);
	}

	bool isIdentity() {
		if (ops().size() == 1 && ops()[0] == "I") {
			return true;
		} else {
			return false;
		}
	}

	std::complex<double>& coeff() {
		return std::get<0>(*this);
	}

	Term& operator*=( const Term& v ) noexcept;

};

class PauliOperator: public tao::operators::commutative_ring<PauliOperator>,
		public tao::operators::equality_comparable<PauliOperator>,
		public tao::operators::commutative_multipliable<PauliOperator, double>,
		public tao::operators::commutative_multipliable<PauliOperator,
				std::complex<double>> {
protected:

	std::unordered_map<std::string, Term> terms;

public:

	PauliOperator();
	PauliOperator(std::complex<double> c);
	PauliOperator(double c);
	PauliOperator(const PauliOperator& i);
	PauliOperator(std::map<int, std::string> operators);
	PauliOperator(std::map<int, std::string> operators, std::string var);
	PauliOperator(std::map<int, std::string> operators,
			std::complex<double> coeff);
	PauliOperator(std::map<int, std::string> operators,
			double coeff);
	PauliOperator(std::map<int, std::string> operators,
			std::complex<double> coeff, std::string var);
	PauliOperator(std::vector<std::map<int, std::string>> operators);

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

	PauliOperator& operator+=( const PauliOperator& v ) noexcept;
	PauliOperator& operator-=( const PauliOperator& v ) noexcept;
	PauliOperator& operator*=( const PauliOperator& v ) noexcept;
	bool operator==( const PauliOperator& v ) noexcept;
	PauliOperator& operator*=( const double v ) noexcept;
	PauliOperator& operator*=( const std::complex<double> v ) noexcept;
};
}
}

bool operator==(const xacc::vqe::PauliOperator& lhs,
	const xacc::vqe::PauliOperator& rhs);

#endif
