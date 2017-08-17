
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
#define BOOST_TEST_MODULE JordanWignerIRTransformationTester

#include <boost/test/included/unit_test.hpp>
#include "JordanWignerIRTransformation.hpp"

using namespace xacc::vqe;

BOOST_AUTO_TEST_CASE(checkJWTransform) {

	auto Instruction = std::make_shared<FermionInstruction>(
			std::vector<std::pair<int, int>> { { 2, 1 }, { 0, 0 }}, 3.17);
	auto Instruction2 = std::make_shared<FermionInstruction>(
			std::vector<std::pair<int, int>> { { 0, 1 }, { 2, 0 }}, 3.17);

	auto kernel = std::make_shared<FermionKernel>("foo");
	kernel->addInstruction(Instruction);
	kernel->addInstruction(Instruction2);

	auto ir = std::make_shared<FermionIR>();
	ir->addKernel(kernel);

	JordanWignerIRTransformation jwTransform;

	auto newIr = jwTransform.transform(ir);

	BOOST_VERIFY(
			"(1.585,0) * X0 * Z1 * X2 + (1.585,0) * Y0 * Z1 * Y2"
					== jwTransform.getResult().toString(""));

	for (auto k : newIr->getKernels()) {
		for (auto i : k->getInstructions()) {
			std::cout << "InstTest: " << i->getName() << "\n";
		}
	}
}

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

	// First off, split the string into lines
	std::vector<std::string> lines, fLineSpaces;
	boost::split(lines, code, boost::is_any_of("\n"));
	auto functionLine = lines[0];
	std::cout << "HELLO WORLD: " << functionLine << "\n";
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


	JordanWignerIRTransformation t;
	auto ir = t.transform(fermionir);

	auto result = t.getResult();
	auto resultsStr = result.toString("");
	boost::replace_all(resultsStr, "+", "+\n");

	std::string expected = R"expected((-0.490661,0) * I +
 (0.0939372,0) * Z1 +
 (0.0939372,0) * Z0 +
 (0.138381,0) * Z0 * Z1 +
 (-0.0365278,0) * Z2 +
 (0.082831,0) * Z0 * Z2 +
 (-0.0365278,0) * Z3 +
 (0.082831,0) * Z1 * Z3 +
 (0.146066,0) * Z2 * Z3 +
 (0.1401,0) * Z1 * Z2 +
 (-0.0572692,0) * X0 * X1 * Y2 * Y3 +
 (0.0572692,0) * X0 * Y1 * Y2 * X3 +
 (0.0572692,0) * Y0 * X1 * X2 * Y3 +
 (-0.0572692,0) * Y0 * Y1 * X2 * X3 +
 (0.1401,0) * Z0 * Z3)expected";

	std::cout << "HELLO: " << resultsStr << "\n";
	BOOST_VERIFY(resultsStr == expected);

	int i = 0;
	for (auto k : ir->getKernels()) {
		std::cout << "KERNEL: " << i << "\n";
		for (auto i : k->getInstructions()) {
			std::cout << "\tInst: " << i->getName() << "\n";
		}
		i++;
	}

}
