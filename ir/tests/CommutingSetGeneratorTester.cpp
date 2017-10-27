
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

	auto inst1 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "Y" },
			{ 1, "Z" }, { 2, "X" } });
	auto inst2 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "X" },
			 { 1, "Z" }, {2, "Y"} });
	auto inst3 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 1, "Y" },
			{ 2, "Z" }, { 3, "X"} });
	auto inst4 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 1, "X" },
			{ 2, "Z" }, {3, "Y"} });
	auto inst5 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "Y" },
			{ 1, "Y" }, {2, "Y"}, {3, "X"} });
	auto inst6 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "Y" },
			{ 1, "X" }, {2, "Y"}, {3, "Y"} });
	auto inst7 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "X" },
			{ 1, "X" }, {2, "Y"}, {3, "X"} });
	auto inst8 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "X" },
			{ 1, "Y" }, {2, "Y"}, {3, "Y"} });
	auto inst9 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "Y" },
			{ 1, "X" }, {2, "X"}, {3, "X"} });
	auto inst10 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "Y" },
				{ 1, "Y" }, {2, "X"}, {3, "Y"} });
	auto inst11 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "X" },
				{ 1, "Y" }, {2, "X"}, {3, "X"} });
	auto inst12 = std::make_shared<SpinInstruction>(std::vector<std::pair<int, std::string>> { { 0, "X" },
				{ 1, "X" }, {2, "X"}, {3, "Y"} });

	CompositeSpinInstruction composite;
	composite.addInstruction(inst3);
	composite.addInstruction(inst11);
	composite.addInstruction(inst9);
	composite.addInstruction(inst1);
	composite.addInstruction(inst2);
	composite.addInstruction(inst5);
	composite.addInstruction(inst12);
	composite.addInstruction(inst10);
	composite.addInstruction(inst6);
	composite.addInstruction(inst8);
	composite.addInstruction(inst7);
	composite.addInstruction(inst4);

	auto str = composite.toString("");
	boost::replace_all(str, "+", "+\n");
	std::cout << "OP:\n" << str << "\n";
	CommutingSetGenerator gen;
	auto sets = gen.getCommutingSet(composite, 4);

	BOOST_VERIFY(sets.size() == 2);
	BOOST_VERIFY(sets[0][0] == 0);
	BOOST_VERIFY(sets[0][1] == 2);
	BOOST_VERIFY(sets[0][2] == 7);
	BOOST_VERIFY(sets[0][3] == 9);
	BOOST_VERIFY(sets[0][4] == 10);
	BOOST_VERIFY(sets[0][5] == 11);

	BOOST_VERIFY(sets[1][0] == 1);
	BOOST_VERIFY(sets[1][1] == 3);
	BOOST_VERIFY(sets[1][2] == 4);
	BOOST_VERIFY(sets[1][3] == 5);
	BOOST_VERIFY(sets[1][4] == 6);
	BOOST_VERIFY(sets[1][5] == 8);


}

