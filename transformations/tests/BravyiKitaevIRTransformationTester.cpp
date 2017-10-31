
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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BravyiKitaevIRTransformationTester

#include <boost/test/included/unit_test.hpp>
#include "BravyiKitaevIRTransformation.hpp"
#include "XACC.hpp"

using namespace xacc::vqe;


std::shared_ptr<FermionKernel> compileKernel(const std::string src) {
	// First off, split the string into lines
	std::vector<std::string> lines, fLineSpaces;
	boost::split(lines, src, boost::is_any_of("\n"));
	auto functionLine = lines[0];
//	std::cout << "HELLO WORLD: " << functionLine << "\n";
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

	xacc::setOption("n-qubits", "4");
	return fermionKernel;
}

/*BOOST_AUTO_TEST_CASE(checkBKTransform) {

	xacc::setOption("n-qubits", "3");
	auto Instruction = std::make_shared<FermionInstruction>(
			std::vector<std::pair<int, int>> { { 2, 1 }, { 0, 0 }}, 3.17);
	auto Instruction2 = std::make_shared<FermionInstruction>(
			std::vector<std::pair<int, int>> { { 0, 1 }, { 2, 0 }}, 3.17);

	auto kernel = std::make_shared<FermionKernel>("foo");
	kernel->addInstruction(Instruction);
	kernel->addInstruction(Instruction2);

	auto ir = std::make_shared<FermionIR>();
	ir->addKernel(kernel);

	BravyiKitaevIRTransformation bkTransform;

	auto newIr = bkTransform.transform(ir);

	auto result = bkTransform.getResult();
	std::cout << "HI: " << bkTransform.getResult().toString("") << "\n";
	BOOST_VERIFY(
			"(-1.585,0) * Y0 * Y1 + (-1.585,0) * X0 * X1 * Z2"
					== bkTransform.getResult().toString(""));

	for (auto k : newIr->getKernels()) {
		for (auto i : k->getInstructions()) {
			std::cout << "InstTest: " << i->getName() << "\n";
		}
	}
}*/

BOOST_AUTO_TEST_CASE(checkH2Transform) {

	const std::string code = R"code(__qpu__ H2_sto-3g_singlet_H2_Molecule() {
        0.276761058984 1 1 0 1 0 0 1 0
        0.280200438869 2 1 0 1 0 0 2 0
        -0.910498000093 1 1 1 0
        0.280200438869 3 1 1 1 1 0 3 0
        0.292131000327 3 1 2 1 2 0 3 0
        0.280200438869 1 1 2 1 2 0 1 0
        -0.664937888133 2 1 2 0
        0.114538342894 2 1 1 1 3 0 0 0
        0.280200438869 1 1 3 1 3 0 1 0
        0.276761058984 0 1 1 1 1 0 0 0
        0.114538342894 3 1 2 1 0 0 1 0
        0.292131000327 2 1 3 1 3 0 2 0
        0.114538342894 2 1 0 1 2 0 0 0
        0.114538342894 3 1 0 1 2 0 1 0
        -0.910498000093 0 1 0 0
        0.114538342894 0 1 2 1 0 0 2 0
        0.114538342894 3 1 1 1 3 0 1 0
        -0.664937888133 3 1 3 0
        0.280200438869 3 1 0 1 0 0 3 0
        0.280200438869 0 1 2 1 2 0 0 0
        0.114538342894 0 1 3 1 1 0 2 0
        0.114538342894 1 1 2 1 0 0 3 0
        0.280200438869 2 1 1 1 1 0 2 0
        0.114538342894 0 1 1 1 3 0 2 0
        0.114538342894 1 1 3 1 1 0 3 0
        0.354466478596
        0.114538342894 2 1 3 1 1 0 0 0
        0.114538342894 1 1 0 1 2 0 3 0
        0.280200438869 0 1 3 1 3 0 0 0
})code";

	auto fermionKernel = compileKernel(code);

	// Create the FermionIR to pass to our transformation.
	auto fermionir = std::make_shared<FermionIR>();
	fermionir->addKernel(fermionKernel);

	xacc::setOption("n-qubits", "4");

	BravyiKitaevIRTransformation t;
	auto ir = t.transform(fermionir);

	auto result = t.getResult();

	BOOST_VERIFY(result.nInstructions() == 15);

	std::vector<std::vector<std::pair<int, std::string>>> expectedTerms {
		{{0,"Z"}},
		{{1,"Z"}},
		{{2,"Z"}},
		{{0,"Z"},{1,"Z"}},
		{{0,"Z"},{2,"Z"}},
		{{1,"Z"},{3,"Z"}},
		{{0,"X"},{1,"Z"},{2,"X"}},
		{{0,"Y"},{1,"Z"},{2,"Y"}},
		{{0,"Z"},{1,"Z"},{2,"Z"}},
		{{0,"Z"},{2,"Z"},{3,"Z"}},
		{{1,"Z"},{2,"Z"},{3,"Z"}},
		{{0,"X"},{1,"Z"},{2,"X"},{3,"Z"}},
		{{0,"Y"},{1,"Z"},{2,"Y"},{3,"Z"}},
		{{0,"Z"},{1,"Z"},{2,"Z"},{3,"Z"}},
		{{0,"I"}}
	};

	for (auto inst : result.getInstructions()) {
		auto cast = std::dynamic_pointer_cast<SpinInstruction>(inst);
		auto terms = cast->getTerms();
		BOOST_VERIFY(std::find(expectedTerms.begin(), expectedTerms.end(), terms) != expectedTerms.end());
	}


}

