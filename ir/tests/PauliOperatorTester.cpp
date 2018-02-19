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
#define BOOST_TEST_MODULE PauliOperatorTester

#include <boost/test/included/unit_test.hpp>
#include "PauliOperator.hpp"
#include <boost/algorithm/string.hpp>
#include "XACC.hpp"
#include "GateInstruction.hpp"
#include "GateFunction.hpp"
#include "GateQIR.hpp"

using namespace xacc::vqe;

BOOST_AUTO_TEST_CASE(checkEasy) {

//	CURRENT: (0.25,0) X1 Z2 X3 + (0,0.25) X1 Z2 Y3 + (0.25,0) Y1 Z2 Y3 + (0,-0.25) Y1 Z2 X3
//	SUM: (0,-0.5) Z0 Y1 + (0.5,0) Z0 X1

	PauliOperator op({{1, "X"}, {2, "Z"}, {3,"X"}}, .25);
	op += PauliOperator({{1,"X"}, {2, "Z"}, {3,"Y"}}, std::complex<double>(0,.25));
	op += PauliOperator({{1,"Y"}, {2, "Z"}, {3,"Y"}}, .25);
	op += PauliOperator({{1,"Y"}, {2, "Z"}, {3,"X"}}, std::complex<double>(0,-.25));

	PauliOperator sum({{0,"Z"}, {1, "Y"}}, std::complex<double>(0,-.5));
	sum += PauliOperator({{0, "Z"}, {1,"X"}}, .5);

	std::cout << "PRINT: " << op.toString() << "\n";
	std::cout << "SUM: " << sum.toString() << "\n";

	op *= sum;

	std::cout << "CURRENTMULT: " << op.toString() << "\n";
}
BOOST_AUTO_TEST_CASE(checkConstruction) {

	PauliOperator inst({ { 4, "X" },
			{ 3, "Z" }, { 9, "Y" }, { 1, "Z" } });

	std::cout << "PRINT: " << inst.toString() << "\n";

	PauliOperator i2({ { 4, "X" },
			{ 3, "Z" } });
	auto sumInst = inst + i2;
	std::cout << "PRINT: " << sumInst.toString() << "\n";

	PauliOperator expected({{3,"Z"}, {4,"X"}});
	expected += PauliOperator({{1,"Z"}, {3,"Z"},{4,"X"}, {9,"Y"}});
	std::cout << "EXPECTED: " << expected.toString() << "\n";
	BOOST_VERIFY(expected == sumInst);
	BOOST_VERIFY(sumInst == sumInst);

	PauliOperator i3({ { 3, "X" } });
	auto testPauliProducts = i2 * i3;
	std::cout << testPauliProducts.toString() << "\n";

	expected.clear();
	expected += PauliOperator({{3,"Y"},{4,"X"}}, std::complex<double>(0,1));
//	expected.addTerm({{3,"Y"},{4,"X"}}, std::complex<double>(0,1));

	std::cout << "NEWEXPECTED: " << expected.toString() << "\n";
	BOOST_VERIFY(expected == testPauliProducts);

	PauliOperator testEqual({ { 4, "X" },
			{ 3, "Z" } }, std::complex<double>(33.3,0.0));
	BOOST_VERIFY(i2 == testEqual);
	BOOST_VERIFY(i2 != inst);

	// Term * scalar multiple
	auto multScalar = testPauliProducts * 3.3;
	BOOST_VERIFY("(0,3.3) Y3 X4" == multScalar.toString());
	BOOST_VERIFY(PauliOperator( { {3,"Y"}, {4,"X"}}, std::complex<double>(0,3.3)) == multScalar);

	testPauliProducts*=3.3;
	BOOST_VERIFY(testPauliProducts == multScalar);

	// Plus with same terms
	sumInst = inst + 2.2 * inst;
	std::cout << inst.toString() << ", " << (2.2*inst).toString() << ", " << sumInst.toString() << "\n";
	BOOST_VERIFY(PauliOperator( { {1,"Z"}, {3,"Z"}, {4,"X"}, {9,"Y"}}, std::complex<double>(3.2,0)) == sumInst);


	PauliOperator i4( { { 0, "Z" } });
	PauliOperator i6( { { 1, "Y" } }, std::complex<double>(0,1));
	i6 += PauliOperator({{1,"X"}});
	PauliOperator i7( { { 1, "Z" } });

	expected.clear();
	expected += PauliOperator({{0, "Z"}, {1, "Y"}}, std::complex<double>(0,1));
	expected += PauliOperator({{0, "Z"}, {1, "X"}});
//	expected.addTerm({{0, "Z"}, {1, "Y"}}, std::complex<double>(0,1));
//	expected.addTerm({{0, "Z"}, {1, "X"}});

	std::cout << i4.toString() << ", " << i6.toString() << "\nSTARTING:\n";
	auto newMultByComp = i4 * i6;
	std::cout << "ENDED:\nHI: " << newMultByComp.toString() << "\n";
	BOOST_VERIFY(expected == newMultByComp);

	std::cout << "MADE IT HERE\n";
	// Z1 * X1 + i * Z1 * Y1 = i * Y1 + i * (-i * X1)
	newMultByComp = i7 * i6;
	expected.clear();
	expected += PauliOperator({{1,"X"}}) + PauliOperator({{1,"Y"}},std::complex<double>(0,1) );
//	expected.addTerm({{1,"X"}}).addTerm({{1,"Y"}}, std::complex<double>(0,1));

	BOOST_VERIFY(expected == newMultByComp);

}

BOOST_AUTO_TEST_CASE(checkVariableCoefficient) {
	// make a4dag a3dag a9 a1
	PauliOperator inst({ { 4, "X" },
			{ 3, "Z" }, }, "theta");
	auto three = 3*inst;
	auto sum = ((inst+inst) + three);
	std::cout << "SUMSTR: " << sum.toString() << "\n";
	BOOST_VERIFY(sum.toString() == "(5,0) theta Z3 X4");

	PauliOperator i2({ { 4, "X" },
			{ 3, "Z" } }, "theta2");

	BOOST_VERIFY(i2 != inst);
	auto sumInst = inst + i2;
	std::cout << "SUMSTR: " << sumInst.toString() << "\n";
	BOOST_VERIFY("(1,0) theta2 Z3 X4 + (1,0) theta Z3 X4" == sumInst.toString());
}

BOOST_AUTO_TEST_CASE(checkBinaryVector) {
//	SpinInstruction inst( { { 4, "X" },
//				{ 3, "Z" }, { 9, "Y" }, { 1, "Z" } });
//
//	// This is Z1 Z3 X4 Y9
//
//	auto bv = inst.toBinaryVector(10);
//
//	BOOST_VERIFY(bv.size() == 20);
//
//	std::vector<int> expected(20);
//	expected[4] = 1;
//	expected[9] = 1;
//	expected[19] = 1;
//	expected[13] = 1;
//	expected[11] = 1;
//
//	BOOST_VERIFY(expected == bv);
//
//	SpinInstruction i;
//	std::vector<int> zeros(20);
//	BOOST_VERIFY(zeros == i.toBinaryVector(10));
//
//	i.fromBinaryVector(expected, std::complex<double>(1.0,0.0));
//
//	BOOST_VERIFY(inst.toString() == i.toString());
}

BOOST_AUTO_TEST_CASE(checkAction) {

	PauliOperator inst({ { 0, "X" },
					{ 1, "Z" }, { 2, "Y" }, { 3, "Z" } }, std::complex<double>(2.2,0));

	auto result = inst.computeActionOnKet("0110");

	BOOST_VERIFY(result[0].first == "1100");
	BOOST_VERIFY(result[0].second == std::complex<double>(0,2.2));

	std::cout << "HELLO:\n" << result[0].first << ", " << result[0].second << "\n";

	PauliOperator inst2({ { 0, "X" },
					{ 1, "Z" }, { 2, "Y" }, { 3, "Y" } }, std::complex<double>(2.2,0));

	auto result2 = inst2.computeActionOnKet("0110");

	std::cout << "HELLO:\n" << result2[0].first << ", " << result2[0].second << "\n";

	BOOST_VERIFY(result2[0].first == "1101");
	BOOST_VERIFY(result2[0].second == std::complex<double>(-2.2,0));

}

BOOST_AUTO_TEST_CASE(checkComplexTerms) {
//	SpinInstruction inst({{1,"X"}, {3,"Y"}, {8,"Z"}}, std::complex<double>(0.5,0.0)),
//			inst2({{1,"Z"}, {3,"X"}, {8,"Z"}}, std::complex<double>(1.2,0.0)),
//			inst3({{1,"Z"}, {3,"Y"}, {9,"Z"}}, std::complex<double>(0.0,1.4)),
//			c({}, std::complex<double>(.25+1.2*1.2, 0.0)
//					+ std::complex<double>(0,1.4)*std::complex<double>(0,1.4)),
//			c2({{1,"Y"},{3,"Z"}}, std::complex<double>(0,2)*std::complex<double>(0,1)*std::complex<double>(.5,0)*std::complex<double>(1.2,0));
//
//	auto op = inst+inst2+inst3;
//
//	auto res = op*op;
//	std::cout << "RES: " << res.toString("") << "\n";
//	res.simplify();
//	std::cout << "Simplified: " << res.toString("") << "\n";
//
//	auto t = c + c2;
//	t.simplify();
//	std::cout << "TSimp: " << t.toString("") << "\n";
//
//	BOOST_VERIFY(t == res);

}


BOOST_AUTO_TEST_CASE(checkHelpMe) {

	PauliOperator sum1({ { 0, "X" } },
				std::complex<double>(.5, 0));
	PauliOperator sum2({ { 0, "Y" } },
					std::complex<double>(0, 0.5));
	PauliOperator current1({ { 0, "X" } },
				std::complex<double>(-.62215, 0));
	PauliOperator current2({ { 0, "Y" } },
					std::complex<double>(0, 0.62215));

	PauliOperator sum, current;
	sum += sum1 + sum2;
	current += current1 + current2;


//	UM: (0.5,0) * X0 + (0,0.5) * Y0
//	Current: (-0.62215,0) * X0 + (0,0.62215) * Y0

	std::cout << "SUM: " << sum.toString() << "\n";
	std::cout << "Current: " << current.toString() << "\n";
	auto x = current * sum;

	std::cout << "RESULT: " << x.toString() << "\n";

	std::cout << (std::complex<double>(1,0) * std::complex<double>(0,1)) << "\n";
}

BOOST_AUTO_TEST_CASE(checkIdentity) {
	PauliOperator i1(
			std::complex<double>(3.3, 0));
	PauliOperator i2(
			std::complex<double>(5.3, 0));
	PauliOperator compInst;
	compInst += i1;

	compInst = compInst + i2;

	std::cout << "HI: " << compInst.toString() << "\n";
	BOOST_VERIFY("(8.6,0)" == compInst.toString());
	BOOST_VERIFY(PauliOperator(std::complex<double>(8.6,0)) == compInst);
}

BOOST_AUTO_TEST_CASE(checkComposition) {
	PauliOperator i1({ { 0, "Z" } });
	PauliOperator i2({ { 1, "X" } });
	PauliOperator i3({ { 1, "Y" } }, std::complex<double>(0,1));
	PauliOperator i4({ { 1, "Z" } });
	PauliOperator i5({ { 1, "Y" } }, std::complex<double>(0,-1));

	PauliOperator compInst, compInst2, compInst3;
	compInst += i2 + i3;

	compInst3 += i2 + i1;

	compInst2 += i2;

	std::complex<double> i(0,1);
	BOOST_VERIFY(compInst == compInst);
	BOOST_VERIFY(compInst != compInst2);
	BOOST_VERIFY(compInst2 == i2);
	BOOST_VERIFY(compInst2 != i1);

	auto multBySpinInst = compInst * i4;
	PauliOperator expected({{1,"X"}}, -1);
	expected += PauliOperator({{1,"Y"}}, -i);
	BOOST_VERIFY(expected == multBySpinInst);

	multBySpinInst = compInst * i1;
	std::cout << "HOWDY: " << multBySpinInst.toString() << "\n";
	expected.clear();
	expected += PauliOperator({{0,"Z"}, {1,"Y"}}, i);
	expected += PauliOperator({{0,"Z"}, {1,"X"}});
	BOOST_VERIFY(expected == multBySpinInst);

	auto multByComposite = compInst * compInst;
	std::cout << "ZERO: " << multByComposite.toString() << "\n";
	BOOST_VERIFY(PauliOperator() == multByComposite);

	auto t = compInst3 * compInst3;
	expected.clear();
	expected += PauliOperator({{0,"Z"}, {1,"X"}}, 2);
	expected += PauliOperator(2);
	BOOST_VERIFY(expected
					== t);

	auto test = compInst * 4.4;
	expected.clear();
	expected += PauliOperator({{1,"X"}}, 4.4);
	expected += PauliOperator({{1,"Y"}}, 4.4*i);
	BOOST_VERIFY(expected == test);

	BOOST_VERIFY(expected == (4.4 * compInst));

	test = compInst * std::complex<double>(4.4, 0);
	BOOST_VERIFY(expected == test);
	BOOST_VERIFY(expected
					== (std::complex<double>(4.4, 0) * compInst));

	// Test Addition now

	auto added = compInst + compInst3;
	expected.clear();
	expected += PauliOperator({{0, "Z"}});
	expected += PauliOperator({{1,"X"}}, 2);
	expected += PauliOperator({{1,"Y"}}, i);

	BOOST_VERIFY(expected == added);
}

BOOST_AUTO_TEST_CASE(checkMatrixElements) {
	PauliOperator op({{0, "X"}, {1, "Y"}, {2, "Z"}});
	auto elements = op.getSparseMatrixElements();
}

BOOST_AUTO_TEST_CASE(checkFromXACCIR) {

	using namespace xacc;
	using namespace xacc::quantum;

	auto f = std::make_shared<GateFunction>("f");
	auto f2 = std::make_shared<GateFunction>("f2");

	auto h1 = GateInstructionRegistry::instance()->create("H", std::vector<int>{0});
	auto h2 = GateInstructionRegistry::instance()->create("H", std::vector<int>{1});
	auto m1 = GateInstructionRegistry::instance()->create("Measure", std::vector<int>{0});
	InstructionParameter p(0), p1(1), p2(2);
	m1->setParameter(0, p);
	auto m2 = GateInstructionRegistry::instance()->create("Measure", std::vector<int>{1});
	m2->setParameter(0,p1);

	auto m3 = GateInstructionRegistry::instance()->create("Measure", std::vector<int>{2});
	m3->setParameter(0,p2);

	f->addInstruction(h1);
	f->addInstruction(h2);
	f->addInstruction(m1);
	f->addInstruction(m2);
	f->addInstruction(m3);

	auto h3 = GateInstructionRegistry::instance()->create("H", std::vector<int>{0});
	auto rx = GateInstructionRegistry::instance()->create("Rx", std::vector<int>{1});
	InstructionParameter q(3.1415/2.0);
	rx->setParameter(0,q);

	auto m4 = GateInstructionRegistry::instance()->create("Measure", std::vector<int>{0});
	InstructionParameter p3(0), p4(1), p5(2);
	m4->setParameter(0, p3);
	auto m5 = GateInstructionRegistry::instance()->create("Measure", std::vector<int>{1});
	m5->setParameter(0,p4);

	auto m6 = GateInstructionRegistry::instance()->create("Measure", std::vector<int>{2});
	m6->setParameter(0,p5);

	f2->addInstruction(h3);
	f2->addInstruction(rx);
	f2->addInstruction(m4);
	f2->addInstruction(m5);
	f2->addInstruction(m6);


	auto ir = std::make_shared<GateQIR>();
	ir->addKernel(f);
	ir->addKernel(f2);

	PauliOperator op;

	op.fromXACCIR(ir);

	std::cout << "HEY: " << op.toString() << "\n";
}

