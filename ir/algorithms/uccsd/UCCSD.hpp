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
#ifndef VQE_IR_ADDUCCSDSTATEPREPARATION_HPP_
#define VQE_IR_ADDUCCSDSTATEPREPARATION_HPP_

#include "IRGenerator.hpp"
#include "FermionKernel.hpp"
#include "FermionIR.hpp"
#include "PauliOperator.hpp"

namespace xacc {

namespace vqe {

/**
 */
class UCCSD: public xacc::IRGenerator {

public:

	/**
	 * Implementations of this method generate a Function IR
	 * instance corresponding to the implementation's modeled
	 * algorithm. The algorithm is specified to operate over the
	 * provided AcceleratorBuffer and can take an optional
	 * vector of InstructionParameters.
	 *
	 * @param bits The bits this algorithm operates on
	 * @return function The algorithm represented as an IR Function
	 */
	virtual std::shared_ptr<Function> generate(
			std::shared_ptr<AcceleratorBuffer> buffer,
			std::vector<InstructionParameter> parameters = std::vector<
					InstructionParameter> { });


	virtual const std::string name() const {
		return "uccsd";
	}

	virtual const std::string description() const {
		return "";
	}
};

}

}

#endif
