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
#ifndef QUANTUM_AQC_FermionKERNEL_HPP_
#define QUANTUM_AQC_FermionKERNEL_HPP_

#include "Function.hpp"
#include "FermionInstruction.hpp"
#include "XACC.hpp"
#include "unsupported/Eigen/CXX11/Tensor"

namespace xacc {
namespace vqe {

/**
 * The FermionKernel is an XACC Function that represents
 * a sum of FermionInstructions.
 */
class FermionKernel: public virtual Function {

protected:

	/**
	 * List of FermionInstructions
	 */
	std::list<InstPtr> instructions;

	/**
	 * This function's name
	 */
	std::string _name;

public:

	/**
	 * The constructor, takes the function unique id and its name.
	 *
	 * @param id
	 * @param name
	 */
	FermionKernel(std::string kernelName) : _name(kernelName) {
	}

	/**
	 * Return the number of FermionInstructions in this sum.
	 *
	 * @return
	 */
	const int nInstructions() override {
		return instructions.size();
	}

    std::shared_ptr<Function> enabledView() override {
        xacc::error("FermionKernel::enabledView is not implemented.");
    }
   void setBits(const std::vector<int> bits)  override {
   }
   void setBitMap(const std::vector<int> bMap) override {}
   const std::vector<int> getBitMap() override {return std::vector<int>{};}

      const int nRequiredBits() const override {return 0;}
   void persist(std::ostream& out) override {};
   void load(std::istream& in) override {};
	/**
	 * Return the FermionInstruction at the given index.
	 *
	 * @param idx The index of the SpinInstruction
	 * @return instruction
	 */
	InstPtr getInstruction(const int idx) override {
		InstPtr ptr;
		if (instructions.size() > idx) {
			ptr = *std::next(instructions.begin(), idx);
		} else {
			xacc::error("Invalid instruction index - " + std::to_string(idx) + ".");
		}
		return ptr;
	}

	void mapBits(std::vector<int> bitMap) override {
	}

	/**
	 * Return all FermionInstruction
	 *
	 * @return instructions
	 */
	std::list<InstPtr> getInstructions() override {
		return instructions;
	}

	/**
	 * Remove the instruction at the given index.
	 *
	 * @param idx The index of the instruction to remove.
	 */
	void removeInstruction(const int idx) override {
		instructions.remove(getInstruction(idx));
	}

	/**
	 * Add an instruction to this quantum
	 * intermediate representation.
	 *
	 * @param instruction
	 */
	void addInstruction(InstPtr instruction) override {
		instructions.push_back(instruction);
	}

	/**
	 * Replace the given current quantum instruction
	 * with the new replacingInst quantum Instruction.
	 *
	 * @param currentInst
	 * @param replacingInst
	 */
	void replaceInstruction(const int idx, InstPtr replacingInst) override {
		std::replace(instructions.begin(), instructions.end(),
				getInstruction(idx), replacingInst);
	}

	/**
	 * Insert the given FermionInstruction at the given index.
	 *
	 * @param idx The index of the new instruction
	 * @param newInst
	 */
	void insertInstruction(const int idx, InstPtr newInst) override {
		auto iter = std::next(instructions.begin(), idx);
		instructions.insert(iter, newInst);
	}

	/**
	 * Return the name of this function
	 * @return
	 */
	const std::string name() const override {
		return _name;
	}

	const std::string description() const override {
		return "";
	}

    std::shared_ptr<Graph> toGraph() override {
        xacc::error("FermionKernel.toGraph not implemented");
    }
    const int depth() override {
        xacc::error("FermionKernel.depth() not implemented.");
    }

    const std::string persistGraph() override {
        xacc::error("FermionKernel.persistGraph() not implemented.");
    }

	/**
	 * Return the qubits this function acts on.
	 * @return
	 */
	const std::vector<int> bits() override {
		return std::vector<int> { };
	}

	const double E_nuc() {
		auto instructions = getInstructions();
		auto instVec = std::vector<InstPtr>(instructions.begin(), instructions.end());

		double e = 0.0;
		// Loop over all Fermionic terms...
		for (int z = 0; z < instructions.size(); ++z) {

			auto f = instVec[z];

			// Get the creation or annihilation sites
			auto termSites = f->bits();
			auto params = f->getParameters();
			auto coeff = params[f->nParameters()-2].as<std::complex<double>>();

			if (termSites.empty()) {
				e = std::real(coeff);
			}
		}
		return e;
	}

	Eigen::Tensor<std::complex<double>, 2> hpq(const int nQubits) {
		Eigen::Tensor<std::complex<double>, 2> hpq(nQubits, nQubits);
		hpq.setZero();

		auto instructions = getInstructions();
		auto instVec = std::vector<InstPtr>(instructions.begin(), instructions.end());

		// Loop over all Fermionic terms...
		for (int z = 0; z < instructions.size(); ++z) {

			auto f = instVec[z];

			// Get the creation or annihilation sites
			auto termSites = f->bits();
			auto params = f->getParameters();
			auto coeff = params[f->nParameters()-2].as<std::complex<double>>();

			if (termSites.size() == 2) {
				int p = termSites[0];
				int q = termSites[1];
				hpq(p,q) = coeff;
			}
		}

		return hpq;
	}

	Eigen::Tensor<std::complex<double>, 4> hpqrs(const int nQubits) {

		Eigen::Tensor<std::complex<double>, 4> hpqrs(nQubits, nQubits, nQubits, nQubits);
		hpqrs.setZero();

		auto instructions = getInstructions();
		auto instVec = std::vector<InstPtr>(instructions.begin(), instructions.end());

		// Loop over all Fermionic terms...
		for (int z = 0; z < instructions.size(); ++z) {

			auto f = instVec[z];

			// Get the creation or annihilation sites
			auto termSites = f->bits();
			auto params = f->getParameters();
			auto coeff = params[f->nParameters()-2].as<std::complex<double>>();

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

	/**
	 * Return an assembly-like string representation for this function .
	 * @param bufferVarName
	 * @return
	 */
	const std::string toString(const std::string& bufferVarName) override {
		std::stringstream ss;
		for (auto i : instructions) {
			ss << i->toString("") << " + \n";
		}
		return ss.str().substr(0, ss.str().size() - 3);
	}

    const std::string toString() override {
        return toString("");
    }

	/**
	 * FermionKernel do not contain parameters. This
	 * method is not implmeneted. An error will be thrown if invoked.
	 *
	 * @param idx
	 * @return
	 */

	InstructionParameter getParameter(const int idx) const override {
		xacc::error("FermionKernel does not contain runtime parameters.");
		return InstructionParameter(0);
	}

	/**
	 * FermionKernel do not contain parameters. This
	 * method is not implmeneted. An error will be thrown if invoked.
	 *
	 * @param idx
	 * @return
	 */

	void setParameter(const int idx, InstructionParameter& p) override {
		xacc::error("FermionKernel does not contain runtime parameters.");
	}

	/**
	 * FermionKernel do not contain parameters. This
	 * method is not implmeneted. An error will be thrown if invoked.
	 *
	 * @param idx
	 * @return
	 */

	std::vector<InstructionParameter> getParameters() override {
		xacc::error("FermionKernel does not contain runtime parameters.");
		return std::vector<InstructionParameter>{};
	}

	void addParameter(InstructionParameter instParam) override {
		xacc::error("FermionKernel does not contain runtime parameters.");
	}

	/**
	 * FermionKernel is not parameterized.
	 */
	bool isParameterized() override {
		return false;
	}

	/**
	 * FermionKernel has 0 parameters.
	 * @return
	 */
	const int nParameters() override {
		return 0;
	}

	/**
	 * FermionKernel do not contain parameters. This
	 * method is not implmeneted. An error will be thrown if invoked.
	 *
	 * @param idx
	 * @return
	 */
	std::shared_ptr<Function> operator()(const std::vector<double>& params) override {
		xacc::error("FermionKernel does not contain runtime parameters.");
		return std::make_shared<FermionKernel>("");
	}
/**
   * Return the number of logical qubits.
   *
   * @return nLogical The number of logical qubits.
   */
  const int nLogicalBits() override {
      return 0;
  }

  /**
   * Return the number of physical qubits.
   *
   * @return nPhysical The number of physical qubits.
   */
  const int nPhysicalBits() override {
      return 0;
  }
  /**
   * Return true if this Instruction has
   * customizable options.
   *
   * @return hasOptions
   */
  bool hasOptions() override {
      return false;
  }

  /**
   * Set the value of an option with the given name.
   *
   * @param optName The name of the option.
   * @param option The value of the option
   */
  void setOption(const std::string optName,
                 InstructionParameter option) override {
      XACCLogger::instance()->error("setOption not implemented for FermionKernel.");
      return;
  }

  /**
   * Get the value of an option with the given name.
   *
   * @param optName Then name of the option.
   * @return option The value of the option.
   */
  InstructionParameter getOption(const std::string optName) override {
       XACCLogger::instance()->error("getOption not implemented for FermionKernel.");
       return InstructionParameter(0);
  }

  /**
   * Return all the Instructions options as a map.
   *
   * @return optMap The options map.
   */
  std::map<std::string, InstructionParameter> getOptions() override {
       XACCLogger::instance()->error("getOptions not implemented for FermionKernel.");
       return std::map<std::string,InstructionParameter>();
  }


	EMPTY_DEFINE_VISITABLE()

};

}
}

#endif
