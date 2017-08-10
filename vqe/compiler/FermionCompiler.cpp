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
#include <regex>
#include <boost/algorithm/string.hpp>
#include "FermionCompiler.hpp"
#include "RuntimeOptions.hpp"
#include "GateQIR.hpp"
#include "ServiceRegistry.hpp"
#include "FermionKernel.hpp"

namespace xacc {

namespace vqe {

std::shared_ptr<IR> FermionCompiler::compile(const std::string& src,
		std::shared_ptr<Accelerator> acc) {
	auto runtimeOptions = RuntimeOptions::instance();

	// Set the Kernel Source code
	kernelSource = src;

	// Here we expect we have a kernel, only one kernel

	// First off, split the string into lines
	std::vector<std::string> lines, fLineSpaces;
	boost::split(lines, src, boost::is_any_of("\n"));
	auto functionLine = lines[0];
	boost::split(fLineSpaces, functionLine, boost::is_any_of(" "));
	auto fName = fLineSpaces[1];
	boost::trim(fName);
	fName = fName.substr(0, fName.find_first_of("("));
	auto firstCodeLine = lines.begin() + 1;
	auto lastCodeLine = lines.end() - 1;
	std::vector<std::string> fermionStrVec(firstCodeLine, lastCodeLine);

	auto fermionKernel = std::make_shared<FermionKernel>("fName");

	int _nQubits = 0;
	// Loop over the lines to create DWQMI
	for (auto termStr : fermionStrVec) {
		boost::trim(termStr);
		if (!termStr.empty() && (std::string::npos != termStr.find_first_of("0123456789"))) {
			std::vector<std::string> splitOnSpaces;
			boost::split(splitOnSpaces, termStr, boost::is_any_of(" "));

			for (auto s : splitOnSpaces) {
				std::cout << s << " ";
			}
			std::cout << "\n";


			// We know first term is coefficient
			// FIXME WHAT IF COMPLEX
			auto coeff = std::stod(splitOnSpaces[0]);
			std::vector<std::pair<int, int>> operators;
			for (int i = 1; i < splitOnSpaces.size()-1; i+=2) {
				auto siteIdx = std::stoi(splitOnSpaces[i]);
				if (siteIdx > _nQubits) {
					_nQubits = siteIdx;
				}
				operators.push_back(
						{siteIdx, std::stoi(
								splitOnSpaces[i + 1]) });
			}

			auto fermionInst = std::make_shared<FermionInstruction>(operators,
					coeff);
			fermionKernel->addInstruction(fermionInst);
		}
	}

	_nQubits++;

	// Create the FermionIR to pass to our transformation.
	auto fermionir = std::make_shared<FermionIR>();
	fermionir->addKernel(fermionKernel);

	// Now we have a Function IR instance that contains information
	// about the fermion representation of the Hamiltonian
	// we are compiling. We need to transform it to a spin
	// hamiltonian.
	std::shared_ptr<IRTransformation> transform;
	if (runtimeOptions->exists("fermion-transformation")) {
		auto transformStr = (*runtimeOptions)["fermion-transformation"];
		transform = ServiceRegistry::instance()->getService<IRTransformation>(
				transformStr);
	} else {
		transform = ServiceRegistry::instance()->getService<IRTransformation>(
				"jordan-wigner");
	}

	// Create the Spin Hamiltonian
	auto transformedIR = transform->transform(fermionir);

	// Prepend State Preparation if requested.
	if (runtimeOptions->exists("state-preparation")) {
		runtimeOptions->insert(std::make_pair("n-qubits", std::to_string(_nQubits)));
		auto statePrepIRTransformStr = (*runtimeOptions)["state-preparation"];
		auto statePrepIRTransform = ServiceRegistry::instance()->getService<
				IRTransformation>(statePrepIRTransformStr);
		return statePrepIRTransform->transform(transformedIR);
	} else {
		return transformedIR;
	}

}

}

}
