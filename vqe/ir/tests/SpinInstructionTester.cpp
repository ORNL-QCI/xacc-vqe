
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
#define BOOST_TEST_MODULE SpinInstructionTester

#include <boost/test/included/unit_test.hpp>
#include "SpinInstruction.hpp"
#include <boost/algorithm/string.hpp>

using namespace xacc::vqe;

BOOST_AUTO_TEST_CASE(checkConstruction) {

	// make a4dag a3dag a9 a1
	SpinInstruction inst(std::vector<std::pair<int, std::string>> { { 4, "X" },
			{ 3, "Z" }, { 9, "Y" }, { 1, "Z" } });

	BOOST_VERIFY(boost::get<std::complex<double>>(inst.coefficient) == std::complex<double>(1,0));
	BOOST_VERIFY(inst.bits().size() == 4);
	BOOST_VERIFY(inst.getParameters().size() == 1);

	BOOST_VERIFY((inst.bits() == std::vector<int> {1, 3, 4, 9}));
	BOOST_VERIFY((
			inst.getParameters() == std::vector<xacc::InstructionParameter> {
						xacc::InstructionParameter(std::complex<double>(1.0, 0.0))}));

	SpinInstruction i2(std::vector<std::pair<int, std::string>> { { 4, "X" },
			{ 3, "Z" } });
	auto sumInst = inst + i2;
	std::string expected = "(1,0) * Z1 * Z3 * X4 * Y9 + (1,0) * Z3 * X4";
	BOOST_VERIFY(expected == sumInst.toString(""));

	SpinInstruction i3(std::vector<std::pair<int, std::string>> { { 3, "X" } });
	auto testPauliProducts = i2 * i3;
	BOOST_VERIFY("(0,1) * Y3 * X4" == testPauliProducts.toString(""));

	BOOST_VERIFY(inst == inst);

	SpinInstruction testEqual(std::vector<std::pair<int, std::string>> { { 4, "X" },
			{ 3, "Z" } }, std::complex<double>(33.3,0.0));

	BOOST_VERIFY(i2 == testEqual);
	BOOST_VERIFY(i2 != inst);

	// Term * scalar multiple
	auto multScalar = testPauliProducts * 3.3;
	BOOST_VERIFY("(0,3.3) * Y3 * X4" == multScalar.toString(""));

	testPauliProducts*=3.3;
	BOOST_VERIFY(testPauliProducts == multScalar);

	// Plus with same terms
	sumInst = inst + 2.2 * inst;
	BOOST_VERIFY("(3.2,0) * Z1 * Z3 * X4 * Y9" == sumInst.toString(""));

	SpinInstruction i4(std::vector<std::pair<int, std::string>> { { 0, "Z" } });
	SpinInstruction i5(std::vector<std::pair<int, std::string>> { { 1, "X" } });
	SpinInstruction i6(std::vector<std::pair<int, std::string>> { { 1, "Y" } }, std::complex<double>(0,1));
	SpinInstruction i7(std::vector<std::pair<int, std::string>> { { 1, "Z" } });

	CompositeSpinInstruction compInst;
	compInst.addInstruction(std::make_shared<SpinInstruction>(i5));
	compInst.addInstruction(std::make_shared<SpinInstruction>(i6));

	auto newMultByComp = i4 * compInst;
	BOOST_VERIFY("(1,0) * Z0 * X1 + (0,1) * Z0 * Y1" == newMultByComp.toString(""));

	// Z1 * X1 + i * Z1 * Y1 = i * Y1 + i * (-i * X1)
	newMultByComp = i7 * compInst;
	BOOST_VERIFY("(0,1) * Y1 + (1,0) * X1" == newMultByComp.toString(""));

}

BOOST_AUTO_TEST_CASE(checkVariableCoefficient) {
	// make a4dag a3dag a9 a1
	SpinInstruction inst(std::vector<std::pair<int, std::string>> { { 4, "X" },
			{ 3, "Z" }, }, "theta");
	auto three = 3*inst;
	auto sum = ((inst+inst) + three);
	BOOST_VERIFY(sum.toString("") == "(5,0) * theta * Z3 * X4");

	SpinInstruction i2(std::vector<std::pair<int, std::string>> { { 4, "X" },
			{ 3, "Z" } }, "theta2");

	BOOST_VERIFY(i2 != inst);
	auto sumInst = inst + i2;
	BOOST_VERIFY("(1,0) * theta * Z3 * X4 + (1,0) * theta2 * Z3 * X4" == sumInst.toString(""));

}
