/***********************************************************************************
 * Copyright (c) 2017, UT-Battelle
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
#include "FCIDumpPreprocessor.hpp"
#include <boost/tokenizer.hpp>
#include "XACC.hpp"

namespace xacc {

namespace vqe {
const std::string FCIDumpPreprocessor::process(const std::string& source,
			std::shared_ptr<Compiler> compiler,
			std::shared_ptr<Accelerator> accelerator) {

	if (boost::contains(source, "FCI") && boost::contains(source, "END")
			&& boost::contains(source, "NELEC")) {

		std::string tempSrc = source;
		std::cout << "Source:\n" << source << "\n";
		auto endIdx = source.find("END");
		auto headerStr = source.substr(0, endIdx + 3);

		std::cout << endIdx << " HEADER:\n" << headerStr << "\n";

		boost::char_separator<char> commasep(",");
		boost::tokenizer<boost::char_separator<char> > commaTokens(headerStr,
				commasep);
		std::vector<std::string> commaSplit;
		std::copy(commaTokens.begin(), commaTokens.end(),
				std::back_inserter<std::vector<std::string> >(commaSplit));

		for (auto s : commaSplit) {
			if (boost::contains(s, "NELEC")) {
				std::vector<std::string> split;
				boost::split(split, s, boost::is_any_of("="));
				boost::trim(split[1]);
				xacc::setOption("n-electrons", split[1]);
			}
		}

		tempSrc.erase(0, endIdx + 4);

		boost::char_separator<char> newlineSep("\n");
		boost::tokenizer<boost::char_separator<char> > linesTokens(tempSrc,
				newlineSep);
		std::vector<std::string> lines;
		std::copy(linesTokens.begin(), linesTokens.end(),
				std::back_inserter<std::vector<std::string> >(lines));

		std::string kernelString = "__qpu__ kernel() {\n";
		for (auto line : lines) {

			std::stringstream buffer(line);
			std::vector<std::string> values {
					std::istream_iterator<std::string>(buffer),
					std::istream_iterator<std::string>() };

			double coeff = std::stod(values[0]);
			int i = std::stoi(values[1]);
			int a = std::stoi(values[2]);
			int j = std::stoi(values[3]);
			int b = std::stoi(values[4]);

			if (i == 0 && j == 0 && a == 0 && b == 0) {
				kernelString += "   " + values[0] + "\n";
			} else if (a == 0 && j == 0 && b == 0) {
				// DO NOTHING FOR NOW
			} else if (j == 0 && b == 0) {
				kernelString += "   " + values[0] + " " + values[1] + " 1 "
						+ values[2] + " 0\n";
			} else {
				kernelString += "   " + values[0] + " " + values[1] + " 1 "
						+ values[3] + " 1 " + values[2] + " 0 " + values[4]
						+ " 0\n";
			}

		}

		kernelString += "}";

		return kernelString;
	} else {
		return source;
	}
}

}

}



