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
#include "AcceleratorDecorator.hpp"
#include "FermionInstruction.hpp"
#include "FermionKernel.hpp"
#include "GateFunction.hpp"
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

TEST(RDMPurificationDecoratorTester, checkGround) {

  if (xacc::hasAccelerator("local-ibm")) {
    xacc::setOption("ibm-shots", "1024");
    xacc::setOption("u-p-depol", ".025");
    xacc::setOption("cx-p-depol", ".03");

    // Get the user-specified Accelerator,
    // or TNQVM if none specified
    auto accelerator = xacc::getAccelerator("local-ibm");
    int nQubits = 4;

    auto accd = xacc::getService<AcceleratorDecorator>("rdm-purification");
    accd->setDecorated(accelerator);

    auto buffer = accd->createBuffer("q", 4);

    xacc::setOption("rdm-source", src);

    // Create the UCCSD ansatz and evaluate
    // at the known optimal angles

    auto compiler = xacc::getService<xacc::Compiler>("xacc-py");

    auto ir = compiler->compile(rucc, accelerator);
    auto ruccsd = ir->getKernel("f");

    Eigen::VectorXd parameters(1);
    parameters << .22984;
    ruccsd = (*ruccsd.get())(parameters);

    auto function = std::make_shared<xacc::quantum::GateFunction>("f");
    function->addInstruction(ruccsd);

    auto buffers =
        accd->execute(buffer, std::vector<std::shared_ptr<Function>>{function});

    std::cout << "Purification Energy: "
              << boost::get<double>(
                     buffers[0]->getInformation("purified-energy"))
              << "\n";

    // EXPECT_NEAR(energy, -1.1371, 1e-4);
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
