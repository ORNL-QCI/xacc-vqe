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
#include "ServiceRegistry.hpp"
#include "FermionKernel.hpp"
#include "MPIProvider.hpp"

namespace xacc {

namespace vqe {

std::shared_ptr<IR> FermionCompiler::compile(const std::string& src,
		std::shared_ptr<Accelerator> acc) {
	auto runtimeOptions = RuntimeOptions::instance();

	// Set the Kernel Source code
	kernelSource = src;

	auto serviceRegistry = xacc::ServiceRegistry::instance();
	std::shared_ptr<MPIProvider> provider;
	if (serviceRegistry->hasService<MPIProvider>("boost-mpi")) {
		provider = serviceRegistry->getService<MPIProvider>("boost-mpi");
	} else {
		provider = serviceRegistry->getService<MPIProvider>("no-mpi");
	}

	provider->initialize();
	auto world = provider->getCommunicator();

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

	fermionKernel = std::make_shared<FermionKernel>("fName");

	for (auto termStr : fermionStrVec) {
		boost::trim(termStr);
		if (!termStr.empty() && (std::string::npos != termStr.find_first_of("0123456789"))) {
			std::vector<std::string> splitOnSpaces;
			boost::split(splitOnSpaces, termStr, boost::is_any_of(" "));

			// We know first term is coefficient
			// FIXME WHAT IF COMPLEX
			auto coeff = std::stod(splitOnSpaces[0]);
			std::vector<std::pair<int, int>> operators;
			for (int i = 1; i < splitOnSpaces.size()-1; i+=2) {
				auto siteIdx = std::stoi(splitOnSpaces[i]);
				if (siteIdx > nQubits) {
					nQubits = siteIdx;
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

	nQubits++;
	xacc::setOption("n-qubits", std::to_string(nQubits));

	// Create the FermionIR to pass to our transformation.
	auto fermionir = std::make_shared<FermionIR>();
	fermionir->addKernel(fermionKernel);

	// Now we have a Function IR instance that contains information
	// about the fermion representation of the Hamiltonian
	// we are compiling. We need to transform it to a spin
	// hamiltonian.
	std::shared_ptr<IRTransformation> transform;
	if (xacc::optionExists("fermion-transformation")) {
		auto transformStr = xacc::getOption("fermion-transformation");
		transform = serviceRegistry->getService<IRTransformation>(
				transformStr);
	} else {
		transform = serviceRegistry->getService<IRTransformation>(
				"jw");
	}

	// Create the Spin Hamiltonian
	if (world->rank() == 0) xacc::info("Mapping Fermion to Spin with " + transform->name());
	auto transformedIR = transform->transform(fermionir);
	if (world->rank() == 0) xacc::info("Done mapping Fermion to Spin.");

	// Prepend State Preparation if requested.
	if (xacc::optionExists("state-preparation")) {
		auto statePrepIRTransformStr = xacc::getOption("state-preparation");
		auto statePrepIRTransform = serviceRegistry->getService<
				IRTransformation>(statePrepIRTransformStr);
		if (world->rank() == 0) xacc::info("Generating State Preparation Circuit with " + statePrepIRTransform->name());
		auto ir = statePrepIRTransform->transform(transformedIR);
		if (world->rank() == 0) xacc::info("Done generating State Preparation Circuit with " + statePrepIRTransform->name());
		return ir;
	} else {
		return transformedIR;
	}

}

}

}
