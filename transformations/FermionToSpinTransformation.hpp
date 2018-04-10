/***********************************************************************************
 * Copyright (c) 2017, UT-Battelle
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
#ifndef VQE_IR_FERMIONTOSPINTRANSFORMATION_HPP_
#define VQE_IR_FERMIONTOSPINTRANSFORMATION_HPP_

#include "IRTransformation.hpp"
#include "FermionKernel.hpp"
#include "FermionIR.hpp"
#include "PauliOperator.hpp"
#include <boost/math/constants/constants.hpp>

#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include "XACC.hpp"
#include "unsupported/Eigen/CXX11/Tensor"

namespace xacc {

namespace vqe {

/**
 */
class FermionToSpinTransformation : public xacc::IRTransformation {

public:

	/**
	 * Return the result of the Jordan Wigner transformation
	 * @return
	 */
	virtual PauliOperator getResult() {
		return result;
	}

	Eigen::Tensor<std::complex<double>, 2> hpq() {
		int nQubits = std::stoi(xacc::getOption("n-qubits"));

		Eigen::Tensor<std::complex<double>, 2> hpq(nQubits, nQubits);
		hpq.setZero();

		auto instructions = fermionKernel->getInstructions();
		auto instVec = std::vector<InstPtr>(instructions.begin(), instructions.end());

		// Loop over all Fermionic terms...
		for (int z = 0; z < instructions.size(); ++z) {

			auto f = instVec[z];

			// Get the creation or annihilation sites
			auto termSites = f->bits();
			auto params = f->getParameters();
			auto coeff = boost::get<std::complex<double>>(params[f->nParameters() - 2]);

			if (termSites.size() == 2) {
				int p = termSites[0];
				int q = termSites[1];
				hpq(p,q) = coeff;
			}
		}

		return hpq;
	}

	Eigen::Tensor<std::complex<double>, 4> hpqrs() {
		int nQubits = std::stoi(xacc::getOption("n-qubits"));

		Eigen::Tensor<std::complex<double>, 4> hpqrs(nQubits, nQubits, nQubits, nQubits);
		hpqrs.setZero();

		auto instructions = fermionKernel->getInstructions();
		auto instVec = std::vector<InstPtr>(instructions.begin(), instructions.end());

		// Loop over all Fermionic terms...
		for (int z = 0; z < instructions.size(); ++z) {

			auto f = instVec[z];

			// Get the creation or annihilation sites
			auto termSites = f->bits();
			auto params = f->getParameters();
			auto coeff = boost::get<std::complex<double>>(params[f->nParameters() - 2]);

			if (termSites.size() == 4) {
				int p = termSites[0];
				int q = termSites[1];
				int r = termSites[2];
				int s = termSites[3];
				hpqrs(p,q,r,s) = coeff;
			}
		}

		return hpqrs;
	}

	virtual PauliOperator transform(FermionKernel& kernel) {
		return PauliOperator();
	}

	bool runParallel = true;

protected:

	/**
	 * Reference to the transformation result.
	 */
	PauliOperator result;

	std::shared_ptr<FermionKernel> fermionKernel;

};

}

}
#endif
