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
#ifndef QUANTUM_AQC_COMPILERS_FermionCOMPILER_HPP_
#define QUANTUM_AQC_COMPILERS_FermionCOMPILER_HPP_

#include "Compiler.hpp"
#include "Utils.hpp"
#include "FermionToSpinTransformation.hpp"
#include "FermionIR.hpp"
#include "unsupported/Eigen/CXX11/Tensor"

namespace xacc {

namespace vqe {

/**
 */
class FermionCompiler: public xacc::Compiler {

public:

	/**
	 * The Compiler.
	 */
	FermionCompiler() {}

	/**
	 *
	 * @param src The QMI source code
	 * @param acc Reference to the D-Wave Accelerator
	 * @return
	 */
	virtual std::shared_ptr<xacc::IR> compile(const std::string& src,
			std::shared_ptr<Accelerator> acc);

	/**
	 *
	 * @return
	 */
	virtual std::shared_ptr<xacc::IR> compile(const std::string& src) {
		return compile(src,nullptr);
	}

	/**
	 * Return the command line options for this compiler
	 *
	 * @return options Description of command line options.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		auto desc = std::make_shared<options_description>(
				"XACC Fermion Compiler Options");
		desc->add_options()("fermion-transformation,T", value<std::string>(),
				"Provide the name of Creation/Annihilation to Spin Transformation. Default is Jordan-Wigner.")(
				"fermion-list-transformations",
				"List all available fermion-to-spin transformations.")
				("no-fermion-transformation", "Skip JW/BK transformation step.")
				("fermion-compiler-silent","Turn off print statements.");
		return desc;
	}

	virtual bool handleOptions(variables_map& args) {
		if (args.count("fermion-list-transformations")) {
			auto ids = xacc::getRegisteredIds<
					IRTransformation>();
			for (auto i : ids) {
				xacc::info("Registered Fermion To Spin Transformation: " + i);
			}
			return true;
		}
		return false;
	}

	/**
	 * We don't allow translations for the DW Compiler.
	 * @param bufferVariable
	 * @param function
	 * @return
	 */
	virtual const std::string translate(const std::string& bufferVariable,
			std::shared_ptr<Function> function) {
		xacc::error("FermionCompiler::translate - Method not implemented");
		return "";
	};

	virtual const std::string name() const {
		return "fermion";
	}

	virtual const std::string description() const {
		return "The Fermion Compiler compiles a high-level second-quantized "
				"fermionic hamiltonian to a spin-hamiltonian amenable for execution "
				"on a QPU.";
	}

	/**
	 * The destructor
	 */
	virtual ~FermionCompiler() {}

protected:

	std::shared_ptr<FermionKernel> fermionKernel;

	int nQubits = 0;

};

}

}

#endif
