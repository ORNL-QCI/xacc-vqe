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
	multBySpinInst.simplify();
//	BOOST_VERIFY("(-1,0) * X1 + (0,-1) * Y1" == multBySpinInst.toString(""));
	multBySpinInst = compInst * i1;
	multBySpinInst.simplify();
	std::cout << "HOWDY: " << multBySpinInst.toString("") << "\n";
	BOOST_VERIFY(
			"(0,1) * Z0 * Y1 + (1,0) * Z0 * X1" == multBySpinInst.toString(""));

	auto multByComposite = compInst * compInst;
	multByComposite.simplify();
	BOOST_VERIFY("(0,0)" == multByComposite.toString(""));

	auto t = compInst3 * compInst3;
	t.simplify();
	BOOST_VERIFY(
			"(2,0) * Z0 * X1 + (2,0) * I"
					== t.toString(""));

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
	added.simplify();
	BOOST_VERIFY("(1,0) * Z0 + (2,0) * X1 + (0,1) * Y1" == added.toString(""));

}

BOOST_AUTO_TEST_CASE(checkIdentity) {
	SpinInstruction i1(std::vector<std::pair<int, std::string>> { { 0, "I" } },
			std::complex<double>(3.3, 0));
	SpinInstruction i2(std::vector<std::pair<int, std::string>> { { 0, "I" } },
			std::complex<double>(5.3, 0));
	CompositeSpinInstruction compInst;
	compInst.addInstruction(std::make_shared<SpinInstruction>(i1));

	compInst = compInst + i2;

	compInst.simplify();

	std::cout << "HI: " << compInst.toString("") << "\n";
	BOOST_VERIFY("(8.6,0) * I" == compInst.toString(""));

}

BOOST_AUTO_TEST_CASE(checkHelpMe) {

	SpinInstruction sum1(std::vector<std::pair<int, std::string>> { { 0, "X" } },
				std::complex<double>(.5, 0));
	SpinInstruction sum2(std::vector<std::pair<int, std::string>> { { 0, "Y" } },
					std::complex<double>(0, 0.5));
	SpinInstruction current1(std::vector<std::pair<int, std::string>> { { 0, "X" } },
				std::complex<double>(-.62215, 0));
	SpinInstruction current2(std::vector<std::pair<int, std::string>> { { 0, "Y" } },
					std::complex<double>(0, 0.62215));

	CompositeSpinInstruction sum, current;
	sum.addInstruction(std::make_shared<SpinInstruction>(sum1));
	sum.addInstruction(std::make_shared<SpinInstruction>(sum2));

	current.addInstruction(std::make_shared<SpinInstruction>(current1));
	current.addInstruction(std::make_shared<SpinInstruction>(current2));

//	UM: (0.5,0) * X0 + (0,0.5) * Y0
//	Current: (-0.62215,0) * X0 + (0,0.62215) * Y0

	std::cout << "SUM: " << sum.toString("") << "\n";
	std::cout << "Current: " << current.toString("") << "\n";
	auto x = current * sum;

	std::cout << "RESULT: " << x.toString("") << "\n";

	std::cout << (std::complex<double>(1,0) * std::complex<double>(0,1)) << "\n";


}
