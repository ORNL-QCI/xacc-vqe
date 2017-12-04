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
#include <boost/mpi.hpp>

#include "chemps2/Hamiltonian.h"
#include "unsupported/Eigen/CXX11/Tensor"

using namespace CheMPS2;

namespace xacc {

namespace vqe {

const std::string FCIDumpPreprocessor::process(const std::string& source,
		std::shared_ptr<Compiler> compiler,
		std::shared_ptr<Accelerator> accelerator) {

	if (boost::contains(source, "FCI") && boost::contains(source, "NELEC")) {
		boost::mpi::communicator world;

		std::string kernelString = "";

		if (world.rank() == 0) {
			int symGroup = 7;
			if (xacc::optionExists("vqe-fcidump-symmetry")) {
				symGroup = std::stoi(xacc::getOption("vqe-fcidump-symmetry"));
				XACCInfo("Setting Symmetry Group to " + std::to_string(symGroup));
			}

			std::ofstream tmpFile(".tmp.fcidump");
			tmpFile << source;
			tmpFile.close();

			Hamiltonian h(".tmp.fcidump", symGroup);

			auto nOrbitals = h.getL();
			auto econst = h.getEconst();

			xacc::setOption("n-qubits", std::to_string(2*nOrbitals));

			Eigen::Tensor<double, 2> hpq(2 * nOrbitals, 2 * nOrbitals);
			Eigen::Tensor<double, 4> hpqrs(2 * nOrbitals, 2 * nOrbitals,
					2 * nOrbitals, 2 * nOrbitals);

			hpq.setZero();
			hpqrs.setZero();

			for (int p = 0; p < nOrbitals; p++) {
				for (int q = 0; q < nOrbitals; q++) {
					hpq(2 * p, 2 * q) = h.getTmat(p, q);
					hpq(2 * p + 1, 2 * q + 1) = h.getTmat(p, q);

					for (int r = 0; r < nOrbitals; r++) {
						for (int s = 0; s < nOrbitals; s++) {

							hpqrs(2 * p, 2 * q + 1, 2 * r + 1, 2 * s) =
									h.getVmat(p, q, s, r) / 2.0;
							hpqrs(2 * p + 1, 2 * q, 2 * r, 2 * s + 1) =
									h.getVmat(p, q, s, r) / 2.0;

							if (p != q && r != s) {
								hpqrs(2 * p, 2 * q, 2 * r, 2 * s) = h.getVmat(p,
										q, s, r) / 2.0;
								hpqrs(2 * p + 1, 2 * q + 1, 2 * r + 1,
										2 * s + 1) = h.getVmat(p, q, s, r)
										/ 2.0;
							}
						}
					}
				}
			}

			std::remove(".tmp.fcidump");

			std::stringstream zz;
			zz << "__qpu__ kernel() {\n" << "   " << std::setprecision(16)
					<< econst << "\n";
			kernelString = zz.str();

			for (int p = 0; p < 2 * nOrbitals; p++) {
				for (int q = 0; q < 2 * nOrbitals; q++) {
					if (std::fabs(hpq(p, q)) > 1e-12) {
						std::stringstream ss;
						ss << std::setprecision(16) << "   " << hpq(p, q) << " "
								<< p << " 1 " << q << " 0\n";
						kernelString += ss.str();
					}

					for (int r = 0; r < 2 * nOrbitals; r++) {
						for (int s = 0; s < 2 * nOrbitals; s++) {
							if (std::fabs(hpqrs(p, q, r, s)) > 1e-12) {
								std::stringstream ss;
								ss << std::setprecision(16) << "   "
										<< hpqrs(p, q, r, s) << " " << p
										<< " 1 " << q << " 1 " << r << " 0 "
										<< s << " 0\n";
								kernelString += ss.str();
							}
						}
					}
				}
			}

			kernelString += "}\n";
		}

		boost::mpi::broadcast(world, kernelString, 0);

		return kernelString;

	} else {

		return source;

	}

}

}

}
