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
#ifndef VQE_IR_COMMONPAULIPRODUCTS_HPP_
#define VQE_IR_COMMONPAULIPRODUCTS_HPP_


#include <map>

namespace xacc {
namespace vqe {

/**
 * CommonPauliProducts keeps track of a map of
 * pairs of pauli products and their equivalent
 * mapping to a coefficient multiplied by a single
 * pauli matrix.
 *
 */
class CommonPauliProducts : public std::map<std::pair<std::string, std::string>,
std::pair<std::complex<double>, std::string>> {

public:

	/**
	 * The constructor
	 */
	CommonPauliProducts() {
		insert(
				std::make_pair(std::make_pair("I", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		insert(
				std::make_pair(std::make_pair("I", "X"),
						std::make_pair(std::complex<double>(1.0, 0.0), "X")));
		insert(
				std::make_pair(std::make_pair("X", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "X")));
		insert(
				std::make_pair(std::make_pair("I", "Y"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Y")));
		insert(
				std::make_pair(std::make_pair("Y", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Y")));
		insert(
				std::make_pair(std::make_pair("Z", "I"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Z")));
		insert(
				std::make_pair(std::make_pair("I", "Z"),
						std::make_pair(std::complex<double>(1.0, 0.0), "Z")));
		insert(
				std::make_pair(std::make_pair("X", "X"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		insert(
				std::make_pair(std::make_pair("Y", "Y"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		insert(
				std::make_pair(std::make_pair("Z", "Z"),
						std::make_pair(std::complex<double>(1.0, 0.0), "I")));
		insert(
				std::make_pair(std::make_pair("X", "Y"),
						std::make_pair(std::complex<double>(0.0, 1.0), "Z")));
		insert(
				std::make_pair(std::make_pair("X", "Z"),
						std::make_pair(std::complex<double>(0.0, -1.0), "Y")));
		insert(
				std::make_pair(std::make_pair("Y", "X"),
						std::make_pair(std::complex<double>(0.0, -1.0), "Z")));
		insert(
				std::make_pair(std::make_pair("Y", "Z"),
						std::make_pair(std::complex<double>(0.0, 1.0), "X")));
		insert(
				std::make_pair(std::make_pair("Z", "X"),
						std::make_pair(std::complex<double>(0.0, 1.0), "Y")));
		insert(
				std::make_pair(std::make_pair("Z", "Y"),
						std::make_pair(std::complex<double>(0.0, -1.0), "X")));
	}
};

}
}

#endif
