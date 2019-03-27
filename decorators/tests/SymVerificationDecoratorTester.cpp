/*******************************************************************************
 * Copyright (c) 2018 UT-Battelle, LLC.
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
#include <gtest/gtest.h>
#include "XACC.hpp"
#include "SymVerificationDecorator.hpp"
#include "PauliOperator.hpp"

using namespace xacc::vqe;

using namespace xacc::quantum;

TEST(SymVerificationDecoratorTester, checkSimple) {
  if (xacc::hasAccelerator("local-ibm")) {
      xacc::setOption("ibm-shots","20000");
    auto acc = xacc::getAccelerator("local-ibm");
    auto buffer = acc->createBuffer("buffer", 2);

    auto compiler = xacc::getService<xacc::Compiler>("xacc-py");
    const std::string src = R"src(def f(buffer, t0):
       Rx(3.1415, 0)
       Ry(1.57079,1)
       Rx(7.8539752,0)
       CNOT(1,0)
       Rz(t0, 0)
       CNOT(1,0)
       Ry(7.8539752, 1)
       Rx(1.57079, 0)
       )src";


    auto ir = compiler->compile(src, acc);
    auto f = ir->getKernel("f");
    f = f->operator()(std::vector<double>(1));

    SymVerificationDecorator decorator;
    decorator.setDecorated(acc);

    PauliOperator op;
    op.fromString("Z0 Z1 + X0 X1 + Y0 Y1 + Z0 + Z1");

    auto measureFunctions = op.toXACCIR()->getKernels();
    for (auto& m : measureFunctions) {
        m->insertInstruction(0,f);
    }

    xacc::setOption("sym-op", "Z0 Z1");

    decorator.execute(buffer, measureFunctions);
  }
}

int main(int argc, char **argv) {
  xacc::Initialize();
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
