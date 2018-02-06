#include "PauliOperator.hpp"
#include "GateQIR.hpp"
#include <boost/math/constants/constants.hpp>

namespace xacc {
namespace vqe {
const std::map<std::string, std::pair<c, std::string>> Term:: pauliProducts =  Term::create_map();

PauliOperator::PauliOperator() {
}

PauliOperator::PauliOperator(std::complex<double> c) {
	terms.emplace(std::make_pair("I",c));
}

PauliOperator::PauliOperator(double c) {
	terms.emplace(std::make_pair("I",c));
}

PauliOperator::PauliOperator(std::string var) {
	terms.emplace(std::make_pair("I", var));
}

PauliOperator::PauliOperator(std::complex<double> c, std::string var) {
	terms.emplace(std::piecewise_construct,
            std::forward_as_tuple("I"),
            std::forward_as_tuple(c, var));
}

PauliOperator::PauliOperator(const PauliOperator& i) : terms(i.terms) {
}

/**
 * The Constructor, takes a vector of
 * qubit-gatename pairs. Initializes coefficient to 1
 *
 * @param operators The pauli operators making up this SpinInstruction
 */
PauliOperator::PauliOperator(std::map<int, std::string> operators) {
	terms.emplace(std::make_pair(Term::id(operators), operators));
}

PauliOperator::PauliOperator(std::map<int, std::string> operators, std::string var) {
	terms.emplace(std::piecewise_construct,
            std::forward_as_tuple(Term::id(operators, var)),
            std::forward_as_tuple(var, operators));
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
	terms.emplace(std::piecewise_construct,
            std::forward_as_tuple(Term::id(operators)),
            std::forward_as_tuple(coeff, operators));

}

PauliOperator::PauliOperator(std::map<int, std::string> operators,
		double coeff) : PauliOperator(operators, std::complex<double>(coeff,0)) {
}

PauliOperator::PauliOperator(std::map<int, std::string> operators,
		std::complex<double> coeff, std::string var) {
	terms.emplace(std::piecewise_construct,
            std::forward_as_tuple(Term::id(operators, var)),
            std::forward_as_tuple(coeff, var, operators));

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

bool PauliOperator::isClose(PauliOperator& other) {
	if (!operator==(other)) {
		return false;
	}

	// These are equal, check their coeffs
	for (auto& kv : terms) {
		Term otherTerm = other.terms[kv.first];
		if (std::abs(kv.second.coeff() - otherTerm.coeff()) > 1e-6) {
			return false;
		}
	}

	return true;
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
		std::complex<double> c = std::get<0>(kv.second);
		std::string v = std::get<1>(kv.second);
		std::map<int, std::string> ops = std::get<2>(kv.second);

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

		auto termId = kv.first;
		auto otherTerm = kv.second;

		if (terms.count(termId)) {
			terms.at(termId).coeff() += otherTerm.coeff();
		} else {
			terms.insert({termId, otherTerm});
		}

		if (std::abs(terms[termId].coeff()) < 1e-12) {
			terms.erase(termId);
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

			if (!newTerms.insert({id, multTerm}).second) {
				newTerms.at(id).coeff() += multTerm.coeff();
			}

			if (std::abs(newTerms.at(id).coeff()) < 1e-12) {
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

	for (auto& kv : terms) {
		bool found = false;
		for (auto& vkv : v.terms) {

			if (kv.second.operator==(vkv.second)) {
				found = true;
				break;
			}
		}

		if (!found) {
			return false;
		}
	}

	return true;
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

	std::string ss;
	std::string myVar = std::get<1>(*this);
	std::string otherVar = std::get<1>(v);
	if (!myVar.empty()) {
		if (!otherVar.empty()) {
			ss = myVar + " " + otherVar;
		} else {
			ss = myVar;
		}
	} else {
		if(!otherVar.empty()) {
			ss = otherVar;
		}
	}

	std::get<1>(*this) = ss;

	auto otherOps = std::get<2>(v);
	for (auto& kv : otherOps) {
		auto qubit = kv.first;
		auto gate = kv.second;
		if (ops().count(qubit)) {
			// This means, we have a op on same qubit in both
			// so we need to check its product
			auto myGate = ops().at(qubit);
			auto gate_coeff = pauliProducts.at(myGate + gate);
			if (gate_coeff.second != "I") {
				ops().at(kv.first) = gate_coeff.second;
			} else {
				ops().erase(kv.first);
			}
			coeff() *= gate_coeff.first;
		} else if (gate != "I") {
			ops().emplace(std::make_pair(qubit, gate));
		}
	}

	return *this;
}

PauliOperator PauliOperator::eval(const std::map<std::string, std::complex<double>> varToValMap) {
	PauliOperator ret;

	for (auto& kv : terms) {

		auto id = kv.first;
		auto term = kv.second;

		for (auto& varVal : varToValMap) {
			if (varVal.first == term.var()) {
				term.var() = "";
				term.coeff() *= varVal.second;
			}
		}

		ret.terms.insert({id, term});
	}

	return ret;
}

std::shared_ptr<IR> PauliOperator::toXACCIR() {
// Create a new GateQIR to hold the spin based terms
	auto newIr = std::make_shared<xacc::quantum::GateQIR>();
	int counter = 0;
	auto pi = boost::math::constants::pi<double>();

	// Populate GateQIR now...
	for (auto& inst : terms) {

		Term spinInst = inst.second;

		// Create a GateFunction and specify that it has
		// a parameter that is the Spin Instruction coefficient
		// that will help us get it to the user for their purposes.

		auto gateFunction = std::make_shared<xacc::quantum::GateFunction>(
				inst.first,
				std::vector<InstructionParameter> { InstructionParameter(
						spinInst.coeff()), InstructionParameter(
						spinInst.isIdentity() ? 1 : 0) });

		// Loop over all terms in the Spin Instruction
		// and create instructions to run on the Gate QPU.
		std::vector<std::shared_ptr<xacc::quantum::GateInstruction>> measurements;
		auto termsMap = spinInst.ops();

		std::vector<std::pair<int, std::string>> terms;
		for (auto& kv : termsMap) {
			if (kv.second != "I" && !kv.second.empty()) {
				terms.push_back( { kv.first, kv.second });
			}
		}

		for (int i = terms.size() - 1; i >= 0; i--) {
			auto qbit = terms[i].first;
			auto gateName = terms[i].second;
			auto gateRegistry =
					xacc::quantum::GateInstructionRegistry::instance();
			auto meas = gateRegistry->create("Measure",
					std::vector<int> { qbit });
			xacc::InstructionParameter classicalIdx(qbit);
			meas->setParameter(0, classicalIdx);
			measurements.push_back(meas);

			if (gateName == "X") {
				auto hadamard = gateRegistry->create("H", std::vector<int> {
						qbit });
				gateFunction->addInstruction(hadamard);
			} else if (gateName == "Y") {
				auto rx = gateRegistry->create("Rx", std::vector<int> { qbit });
				InstructionParameter p(pi / 2.0);
				rx->setParameter(0, p);
				gateFunction->addInstruction(rx);
			}

		}

		if (!spinInst.isIdentity()) {
			for (auto m : measurements) {
				gateFunction->addInstruction(m);
			}
		}

		newIr->addKernel(gateFunction);
		counter++;
	}
	return newIr;
}
}
}

bool operator==(const xacc::vqe::PauliOperator& lhs,
		const xacc::vqe::PauliOperator& rhs) {
	if (lhs.getTerms().size() != rhs.getTerms().size()) {
		return false;
	}
	for (auto& kv : lhs.getTerms()) {
		bool found = false;
		for (auto& vkv : rhs.getTerms()) {

			if (kv.second.operator==(vkv.second)) {
				found = true;
				break;
			}
		}

		if (!found) {
			return false;
		}
	}

	return true;
}
