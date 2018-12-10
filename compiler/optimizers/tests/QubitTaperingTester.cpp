/*******************************************************************************
 * Copyright (c) 2017 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Alexander J. McCaskey - initial API and implementation
 *******************************************************************************/
#include "CNOT.hpp"
#include "Compiler.hpp"
#include "GateFunction.hpp"
#include "GateIR.hpp"
#include "InstructionIterator.hpp"
#include "PauliOperator.hpp"
#include "QubitTapering.hpp"
#include "Rz.hpp"
#include <gtest/gtest.h>

using namespace xacc::vqe;

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

TEST(QubitTaperingTester, checkH2) {

  auto c = xacc::getService<xacc::Compiler>("fermion");
      xacc::info("COMPILING");

  auto ir = c->compile(src);
        xacc::info("COMPILED");

  QubitTapering tapering;
  auto newIR = tapering.transform(ir);

  PauliOperator expected, actual;
  expected.fromString("(-0.335683,0) I + (-0.780643,0) Z0 + (0.181625,0) X0");

  actual.fromXACCIR(newIR);

  EXPECT_TRUE(actual == expected);
}

TEST(QubitTaperingTester, checkH2WithAnsatz) {

  if (xacc::hasAccelerator("local-ibm")) {
    auto acc = xacc::getAccelerator("local-ibm");
    auto b = acc->createBuffer("q",4);
    
    auto compiler = xacc::getService<xacc::Compiler>("xacc-py");
    const std::string src2 = R"src(def f(buffer, t0, t1):
       Ry(t0,0)
       Rz(t1,0)
       )src";

    auto ir = compiler->compile(src2, acc);

    Eigen::VectorXd p(2);
    p(0) = -3.1415;
    p(1) = 1.57;
    auto f = ir->getKernel("f")->operator()(p);

    auto c = xacc::getService<xacc::Compiler>("fermion");
    auto tir = c->compile(src);

    auto measureFunctions = tir->getKernels();
    for (auto &m : measureFunctions) {
      m->insertInstruction(0, f);
    }

    QubitTapering tapering;
    auto newIR = tapering.transform(tir);

    PauliOperator expected, actual;
    expected.fromString("(-0.335683,0) I + (-0.780643,0) Z0 + (0.181625,0) X0");

    actual.fromXACCIR(newIR);

    EXPECT_TRUE(actual == expected);

    std::cout << "HELLO: " << newIR->getKernels()[0]->getInstruction(0)->toString("q") << "\n";

    for (auto& k : newIR->getKernels()) {
        EXPECT_TRUE(k->getInstruction(0)->isComposite());
    }
    
  }
}

int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
