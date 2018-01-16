
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
#define BOOST_TEST_MODULE CommutingSetGeneratorTester

#include <boost/test/included/unit_test.hpp>
#include "CommutingSetGenerator.hpp"

using namespace xacc::vqe;

BOOST_AUTO_TEST_CASE(checkCommutingSets) {

	PauliOperator inst1(std::map<int, std::string> { { 0, "Y" },
			{ 1, "Z" }, { 2, "X" } });
	PauliOperator inst2(std::map<int, std::string> { { 0, "X" },
			 { 1, "Z" }, {2, "Y"} });
	PauliOperator inst3(std::map<int, std::string> { { 1, "Y" },
			{ 2, "Z" }, { 3, "X"} });
	PauliOperator inst4(std::map<int, std::string> { { 1, "X" },
			{ 2, "Z" }, {3, "Y"} });
	PauliOperator inst5(std::map<int, std::string> { { 0, "Y" },
			{ 1, "Y" }, {2, "Y"}, {3, "X"} });
	PauliOperator inst6(std::map<int, std::string> { { 0, "Y" },
			{ 1, "X" }, {2, "Y"}, {3, "Y"} });
	PauliOperator inst7(std::map<int, std::string> { { 0, "X" },
			{ 1, "X" }, {2, "Y"}, {3, "X"} });
	PauliOperator inst8(std::map<int, std::string> { { 0, "X" },
			{ 1, "Y" }, {2, "Y"}, {3, "Y"} });
	PauliOperator inst9(std::map<int, std::string> { { 0, "Y" },
			{ 1, "X" }, {2, "X"}, {3, "X"} });
	PauliOperator inst10(std::map<int, std::string> { { 0, "Y" },
				{ 1, "Y" }, {2, "X"}, {3, "Y"} });
	PauliOperator inst11(std::map<int, std::string> { { 0, "X" },
				{ 1, "Y" }, {2, "X"}, {3, "X"} });
	PauliOperator inst12(std::map<int, std::string> { { 0, "X" },
				{ 1, "X" }, {2, "X"}, {3, "Y"} });

	PauliOperator composite = inst1 + inst2 + inst3 + inst4
			+ inst5 + inst6 + inst7 + inst8 + inst9 + inst10
			+ inst11 + inst12;

	auto str = composite.toString();
	boost::replace_all(str, "+", "+\n");
	std::cout << "OP:\n" << str << "\n";
	CommutingSetGenerator gen;
	auto sets = gen.getCommutingSet(composite, 4);

	std::cout << "SIZE: " << sets.size() << "\n";
//	BOOST_VERIFY(sets.size() == 2);
//	BOOST_VERIFY(sets[0][0] == 0);
//	BOOST_VERIFY(sets[0][1] == 2);
//	BOOST_VERIFY(sets[0][2] == 7);
//	BOOST_VERIFY(sets[0][3] == 9);
//	BOOST_VERIFY(sets[0][4] == 10);
//	BOOST_VERIFY(sets[0][5] == 11);
//
//	BOOST_VERIFY(sets[1][0] == 1);
//	BOOST_VERIFY(sets[1][1] == 3);
//	BOOST_VERIFY(sets[1][2] == 4);
//	BOOST_VERIFY(sets[1][3] == 5);
//	BOOST_VERIFY(sets[1][4] == 6);
//	BOOST_VERIFY(sets[1][5] == 8);


}

