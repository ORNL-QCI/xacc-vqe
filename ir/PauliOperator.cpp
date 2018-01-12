#include "PauliOperator.hpp"

namespace xacc {
namespace vqe {

PauliOperator::PauliOperator() {
}

PauliOperator::PauliOperator(std::complex<double> c) {
	terms["I"] = Term(c);

}
PauliOperator::PauliOperator(double c) {
	terms["I"] = Term(c);
}

/**
 * The copy constructor
 * @param i
 */
PauliOperator::PauliOperator(const PauliOperator& i) : terms(i.terms) {
}

/**
 * The Constructor, takes a vector of
 * qubit-gatename pairs. Initializes coefficient to 1
 *
 * @param operators The pauli operators making up this SpinInstruction
 */
PauliOperator::PauliOperator(std::map<int, std::string> operators) {
	Term t(operators);
	terms[t.id()] = t;
}

PauliOperator::PauliOperator(std::map<int, std::string> operators, std::string var) {
	Term t(var, operators);
	terms[t.id()] = t;
}

/**
 * The Constructor, takes a vector of
 * qubit-gatename pairs and this instruction's coefficient
 *
 * @param operators
 * @param coeff
 */
PauliOperator::PauliOperator(std::map<int, std::string> operators,
		std::complex<double> coeff) {
	Term t(coeff, operators);
	terms[t.id()] = t;
}

PauliOperator::PauliOperator(std::map<int, std::string> operators,
		double coeff) : PauliOperator(operators, std::complex<double>(coeff,0)) {
}

PauliOperator::PauliOperator(std::map<int, std::string> operators,
		std::complex<double> coeff, std::string var) {
	Term t(coeff, var, operators);
	terms[t.id()] = t;
}

PauliOperator::PauliOperator(std::vector<std::map<int, std::string>> operators) {

	for (auto& element : operators) {
		Term t(element);
		auto id = t.id();

		if (terms.count(id)) {
			std::get<0>(terms[id]) += std::complex<double>(1,0);
		} else {
			terms[id] = t;
		}
	}
}

const std::vector<std::pair<std::string, std::complex<double>>> PauliOperator::computeActionOnKet(
		const std::string& bitString) {

	std::vector<std::pair<std::string, std::complex<double>>> ret;
	std::string newBits = bitString;
	std::complex<double> newCoeff(1,0), i(0,1);

	for (auto& kv : terms) {
		newCoeff = kv.second.coeff();
		for (auto& t : std::get<2>(kv.second)) {
			auto idx = t.first;
			auto gate = t.second;
			if (gate == "Z") {
				newCoeff *= newBits[idx] == '1' ? -1 : 1;
			} else if (gate == "X") {
				newBits[idx] = (newBits[idx] == '1' ? '0' : '1');
			} else if (gate == "Y") {
				newCoeff *= newBits[idx] == '1' ? -i : i;
				newBits[idx] = (newBits[idx] == '1' ? '0' : '1');
			}
		}

		ret.push_back(std::make_pair(newBits, newCoeff));
	}

	return ret;
}

const std::vector<std::pair<std::string, std::complex<double>>> PauliOperator::computeActionOnBra(
		const std::string& bitString) {
	std::vector<std::pair<std::string, std::complex<double>>> ret;

	for (auto& kv : terms) {
		std::complex<double> newCoeff(1,0), i(0,1);
		std::string newBits = bitString;
		newCoeff = kv.second.coeff();
		for (auto& t : std::get<2>(kv.second)) {
			auto idx = t.first;
			auto gate = t.second;
			if (gate == "Z") {
				newCoeff *= newBits[idx] == '1' ? -1 : 1;
			} else if (gate == "X") {
				newBits[idx] = (newBits[idx] == '1' ? '0' : '1');
			} else if (gate == "Y") {
				newCoeff *= newBits[idx] == '1' ? i : -i;
				newBits[idx] = (newBits[idx] == '1' ? '0' : '1');
			}
		}

		ret.push_back(std::make_pair(newBits, newCoeff));
	}

	return ret;
}

const int PauliOperator::nTerms() {
	return terms.size();
}

/**
 * Persist this Instruction to an assembly-like
 * string.
 *
 * @param bufferVarName The name of the AcceleratorBuffer
 * @return str The assembly-like string.
 */
const std::string PauliOperator::toString() {
	std::stringstream s;
	for (auto& kv : terms) {
		std::complex<double> c;
		std::string v;
		std::map<int, std::string> ops;
		std::tie(c, v, ops) = kv.second;

		s << c << " ";
		if (!v.empty()) {
			s << v << " ";
		}

		for (auto& kv2 : ops) {
			if (kv2.second == "I") {
				s << "I ";
			} else {
				s << kv2.second << kv2.first << " ";
			}
		}

		s << "+ ";
	}

	auto r = s.str().substr(0, s.str().size() - 2);
	boost::trim(r);
	return r;
}

void PauliOperator::clear() {
	terms.clear();
}

PauliOperator& PauliOperator::operator+=( const PauliOperator& v ) noexcept {
	for (auto& kv : v.terms) {

		if (terms.count(kv.first)) {
			std::get<0>(terms[kv.first]) += std::get<0>(kv.second);
		} else {
			terms[kv.first] = kv.second;
		}

		if (std::abs(terms[kv.first].coeff()) < 1e-12) {
			terms.erase(kv.first);
		}
	}

	return *this;
}

PauliOperator& PauliOperator::operator-=( const PauliOperator& v ) noexcept {
	return operator+=(-1.0 * v);
}

PauliOperator& PauliOperator::operator*=( const PauliOperator& v ) noexcept {

	std::unordered_map<std::string, Term> newTerms;
	for (auto& kv : terms) {
		for (auto& vkv : v.terms) {
			auto multTerm = kv.second * vkv.second;
			auto id = multTerm.id();
			if (newTerms.count(id)) {
				newTerms[id].coeff() += multTerm.coeff();
			} else {
				newTerms[id] = multTerm;
			}

			if (std::abs(newTerms[id].coeff()) < 1e-12) {
				newTerms.erase(id);
			}
		}
	}
	terms = newTerms;
	return *this;
}

bool PauliOperator::operator==( const PauliOperator& v ) noexcept {
	if (terms.size() != v.terms.size()) {
		return false;
	}
	std::vector<std::string> k1, k2;
	for (auto& kv : terms) k1.push_back(kv.first);
	for (auto& kv : v.terms) k2.push_back(kv.first);
	return std::is_permutation(k1.begin(), k1.end(), k2.begin());
}

PauliOperator& PauliOperator::operator*=( const double v ) noexcept {
	return operator*=(std::complex<double>(v,0));
}

PauliOperator& PauliOperator::operator*=( const std::complex<double> v ) noexcept {
	for (auto& kv : terms) {
		std::get<0>(kv.second) *= v;
	}
	return *this;
}

Term& Term::operator*=( const Term& v ) noexcept {

	coeff() *= std::get<0>(v);

	std::stringstream ss;
	if (!std::get<1>(*this).empty()) {
		if (!std::get<1>(v).empty()) {
			ss << std::get<1>(*this) << " " << std::get<1>(v);
		} else {
			ss << std::get<1>(*this);
		}
	} else {
		if(!std::get<1>(v).empty()) {
			ss << std::get<1>(v);
		}
	}

	std::get<1>(*this) = ss.str();

	std::map<int, std::string> newTerms = std::get<2>(*this);
	for (auto& kv : std::get<2>(v)) {
		auto qubit = kv.first;
		auto gate = kv.second;
		if (newTerms.count(qubit)) {
			// This means, we have a op on same qubit in both
			// so we need to check its product
			auto myGate = newTerms[qubit];
			auto gate_coeff = pauliProducts.at(myGate + gate);
			if (gate_coeff.second != "I") {
				newTerms[kv.first] = gate_coeff.second;
			} else {
				newTerms.erase(kv.first);
				if (newTerms.empty()) {
					newTerms[0] = "I";
				}
			}
			coeff() *= gate_coeff.first;
		} else if (gate != "I") {
			newTerms[kv.first] = kv.second;
		}
	}

	std::get<2>(*this) = newTerms;

	return *this;
}
}
}

bool operator==(const xacc::vqe::PauliOperator& lhs,
		const xacc::vqe::PauliOperator& rhs) {
	if (lhs.getTerms().size() != rhs.getTerms().size()) {
		return false;
	}

	std::vector<std::string> k1, k2;
	for (auto& kv : lhs.getTerms())
		k1.push_back(kv.first);
	for (auto& kv : rhs.getTerms())
		k2.push_back(kv.first);

	return std::is_permutation(k1.begin(), k1.end(), k2.begin());
}
