#include "JordanWignerIRTransformation.hpp"
#include "GateQIR.hpp"

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/symbol.h>
#include <symengine/subs.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/parser.h>
#include <symengine/series.h>

using SymEngine::Basic;
using SymEngine::Add;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Symbol;
using SymEngine::symbol;
using SymEngine::umap_basic_num;
using SymEngine::Integer;
using SymEngine::integer;
using SymEngine::multinomial_coefficients;
using SymEngine::RCP;
using SymEngine::rcp_dynamic_cast;
using SymEngine::parse;
using SymEngine::series;
using SymEngine::map_basic_basic;
using SymEngine::subs;

namespace xacc {
namespace vqe {

std::shared_ptr<IR> JordanWignerIRTransformation::transform(
		std::shared_ptr<IR> ir) {

	// We assume we have a FermionIR instance, which contains
	// one FermionKernel, which contains N FermionInstructions, one
	// for each term in the hamiltonian.

	// We want to map that to a Hamiltonian composed of pauli matrices
	// But, we want each term of that to be a separate IR Function.

	// Create a new GateQIR to hold the spin based terms
	auto newIr = std::make_shared<xacc::quantum::GateQIR>();

	auto fermiKernel = ir->getKernels()[0];

	int counter = 1;
	for (auto fermionInst : fermiKernel->getInstructions()) {

		auto gateFunction = std::make_shared<xacc::quantum::GateFunction>(
				"term" + std::to_string(counter));

		auto termSites = fermionInst->bits();
		auto params = fermionInst->getParameters();

		RCP<const Basic> termExpression = SymEngine::real_double(
				boost::get<double>(params[termSites.size()]));

		for (int i = 0; i < termSites.size(); i++) {
			RCP<const Basic> op;
			std::cout << "ADDING SYMBOL " << termSites[i] << ", " << params[i] << "\n";
			// here symbol 'd' is a creation operator
			// and symbol 'a' is an anhililation operator
			if (boost::get<int>(params[i])) { // If this is a creation operator
				op = SymEngine::symbol("d" + std::to_string(termSites[i]));
			} else {
				op = symbol("a" + std::to_string(termSites[i]));
			}

			termExpression = mul(op, termExpression);
		}

		std::cout << "Term: " << *termExpression << "\n";

	}

	std::cout << *parse("3.17 * d0 * a2") << "\n";
	return newIr;
}

}
}

