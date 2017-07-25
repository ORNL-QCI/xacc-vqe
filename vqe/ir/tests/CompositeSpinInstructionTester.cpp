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
#define BOOST_TEST_MODULE CompositeSpinInstructionTester

#include <boost/test/included/unit_test.hpp>
#include "CompositeSpinInstruction.hpp"
#include "SpinInstruction.hpp"

using namespace xacc::vqe;

BOOST_AUTO_TEST_CASE(checkConstruction) {
	SpinInstruction i1(std::vector<std::pair<int, std::string>> { { 0, "Z" } });
	SpinInstruction i2(std::vector<std::pair<int, std::string>> { { 1, "X" } });
	SpinInstruction i3(std::vector<std::pair<int, std::string>> { { 1, "Y" } }, std::complex<double>(0,1));
	SpinInstruction i4(std::vector<std::pair<int, std::string>> { { 1, "Z" } });
	SpinInstruction i5(std::vector<std::pair<int, std::string>> { { 1, "Y" } }, std::complex<double>(0,-1));

	CompositeSpinInstruction compInst, compInst2, compInst3;
	compInst.addInstruction(std::make_shared<SpinInstruction>(i2));
	compInst.addInstruction(std::make_shared<SpinInstruction>(i3));

	compInst3.addInstruction(std::make_shared<SpinInstruction>(i2));
	compInst3.addInstruction(std::make_shared<SpinInstruction>(i1));

	compInst2.addInstruction(std::make_shared<SpinInstruction>(i2));

	BOOST_VERIFY(compInst == compInst);
	BOOST_VERIFY(compInst != compInst2);
	BOOST_VERIFY(compInst2 == i2);
	BOOST_VERIFY(compInst2 != i1);

	auto multBySpinInst = compInst * i4;

	BOOST_VERIFY("(0,-1) * Y1 + (-1,0) * X1" == multBySpinInst.toString(""));
	multBySpinInst = compInst * i1;
	BOOST_VERIFY(
			"(1,0) * X1 * Z0 + (0,1) * Y1 * Z0" == multBySpinInst.toString(""));

	auto multByComposite = compInst * compInst;
	BOOST_VERIFY("(0,0)" == multByComposite.toString(""));

	BOOST_VERIFY(
			"(2,0) * I + (1,0) * X1 * Z0 + (1,0) * Z0 * X1"
					== (compInst3 * compInst3).toString(""));

	auto test = compInst * 4.4;
	BOOST_VERIFY("(4.4,0) * X1 + (0,4.4) * Y1" == test.toString(""));
	BOOST_VERIFY(
			"(4.4,0) * X1 + (0,4.4) * Y1" == (4.4 * compInst).toString(""));

	test = compInst * std::complex<double>(4.4, 0);
	BOOST_VERIFY("(4.4,0) * X1 + (0,4.4) * Y1" == test.toString(""));
	BOOST_VERIFY(
			"(4.4,0) * X1 + (0,4.4) * Y1"
					== (std::complex<double>(4.4, 0) * compInst).toString(""));

	// Test Addition now

	auto added = compInst + compInst3;
	BOOST_VERIFY("(2,0) * X1 + (0,1) * Y1 + (1,0) * Z0" == added.toString(""));

}

