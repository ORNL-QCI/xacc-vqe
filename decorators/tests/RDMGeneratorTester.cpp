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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
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
#include "FermionInstruction.hpp"
#include "FermionKernel.hpp"
#include "IRGenerator.hpp"
#include "Instruction.hpp"
#include "RDMGenerator.hpp"
#include "XACC.hpp"
#include <gtest/gtest.h>

using namespace xacc::vqe;
using namespace xacc;

const std::string src = R"src(__qpu__ kernel() {
   0.7080240949826064
-1.248846801817026 0 1 0 0
-1.248846801817026 1 1 1 0
-0.4796778151607899 2 1 2 0
-0.4796778151607899 3 1 3 0
0.33667197218932576 0 1 1 1 1 0 0 0
0.0908126658307406 0 1 1 1 3 0 2 0
0.09081266583074038 0 1 2 1 0 0 2 0
0.331213646878486 0 1 2 1 2 0 0 0
0.09081266583074038 0 1 3 1 1 0 2 0
0.331213646878486 0 1 3 1 3 0 0 0
0.33667197218932576 1 1 0 1 0 0 1 0
0.0908126658307406 1 1 0 1 2 0 3 0
0.09081266583074038 1 1 2 1 0 0 3 0
0.331213646878486 1 1 2 1 2 0 1 0
0.09081266583074038 1 1 3 1 1 0 3 0
0.331213646878486 1 1 3 1 3 0 1 0
0.331213646878486 2 1 0 1 0 0 2 0
0.09081266583074052 2 1 0 1 2 0 0 0
0.331213646878486 2 1 1 1 1 0 2 0
0.09081266583074052 2 1 1 1 3 0 0 0
0.09081266583074048 2 1 3 1 1 0 0 0
0.34814578469185886 2 1 3 1 3 0 2 0
0.331213646878486 3 1 0 1 0 0 3 0
0.09081266583074052 3 1 0 1 2 0 1 0
0.331213646878486 3 1 1 1 1 0 3 0
0.09081266583074052 3 1 1 1 3 0 1 0
0.09081266583074048 3 1 2 1 0 0 1 0
0.34814578469185886 3 1 2 1 2 0 3 0
})src";

const std::string rucc = R"rucc(def f(buffer, theta):
    X(0)
    X(1)
    Rx(1.5707,0)
    H(1)
    H(2)
    H(3)
    CNOT(0,1)
    CNOT(1,2)
    CNOT(2,3)
    Rz(theta,3)
    CNOT(2,3)
    CNOT(1,2)
    CNOT(0,1)
    Rx(-1.5707,0)
    H(1)
    H(2)
    H(3)
    )rucc";

TEST(RDMGeneratorTester, checkGround) {

  if (xacc::hasAccelerator("tnqvm")) {
    // Get the user-specified Accelerator,
    // or TNQVM if none specified
    auto accelerator = xacc::getAccelerator("tnqvm");
    int nQubits = 4;

    // Compile the above source code to get hpq and hpqrs,
    // and the nuclear energy
    xacc::setOption("fermion-compiler-silent", "");
    xacc::setOption("no-fermion-transformation", "");
    auto c = xacc::getCompiler("fermion");
    auto ir = c->compile(src, accelerator);
    auto fermionKernel =
        std::dynamic_pointer_cast<FermionKernel>(ir->getKernels()[0]);
    xacc::unsetOption("no-fermion-transformation");
    auto energy = fermionKernel->E_nuc();
    auto hpq = fermionKernel->hpq(nQubits);
    auto hpqrs = fermionKernel->hpqrs(nQubits);

    // Create the UCCSD ansatz and evaluate
    // at the known optimal angles
    auto compiler = xacc::getService<xacc::Compiler>("xacc-py");

    auto ir2 = compiler->compile(rucc, accelerator);
    auto ruccsd = ir2->getKernel("f");

    Eigen::VectorXd parameters(1);
    parameters << .22984;
    ruccsd = (*ruccsd.get())(parameters);

    // Create the 2-RDM
    std::vector<int> qubitMap {1,3,5,7}; // map to physical qubits
    ruccsd->mapBits(qubitMap);
    RDMGenerator generator(nQubits, accelerator, hpq, hpqrs);
    auto buffers = generator.generate(ruccsd, qubitMap);

    EXPECT_EQ(buffers.size(), 54);

    // Get the 1 rdm from the 2 rdm
    Eigen::array<int, 2> cc2({1, 3});
    auto rho_pqrs = generator.rho_pqrs;
    Eigen::Tensor<std::complex<double>, 2> rho_pq = rho_pqrs.trace(cc2);

    // Compute the 2 rdm contribution to the energy
    for (int p = 0; p < nQubits; p++) {
      for (int q = 0; q < nQubits; q++) {
        for (int r = 0; r < nQubits; r++) {
          for (int s = 0; s < nQubits; s++) {
            energy +=
                0.5 * std::real(hpqrs(p, q, r, s) *
                                (rho_pqrs(p, q, s, r) + rho_pqrs(r, s, q, p)));
          }
        }
      }
    }
    std::cout << "2rdm e = " << std::setprecision(12) << energy << "\n";

    // Compute the 1 rdm contribution to the energy
    for (int p = 0; p < nQubits; p++) {
      for (int q = 0; q < nQubits; q++) {
        energy += 0.5 * std::real(hpq(p, q) * (rho_pq(p, q) + rho_pq(q, p)));
      }
    }
    std::cout << "ENERGY: " << energy << "\n";

    // auto energy = rhoGen.energy();
    EXPECT_NEAR(energy, -1.1371, 1e-4);
  }
}

int main(int argc, char **argv) {
  xacc::Initialize();
  int ret = 0;
  ::testing::InitGoogleTest(&argc, argv);
  ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
