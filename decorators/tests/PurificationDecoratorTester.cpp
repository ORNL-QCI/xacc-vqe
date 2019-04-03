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
#include "PurificationDecorator.hpp"
#include "PauliOperator.hpp"
#include "xacc_service.hpp"

using namespace xacc;
using namespace xacc::quantum;

using namespace xacc::vqe;

TEST(ImprovedSamplingDecoratorTester, checkMultiple) {
  int shots = 8192;
  int nExecs = 2;

  if (xacc::hasAccelerator("local-ibm")) {
    auto acc = xacc::getAccelerator("local-ibm");
    xacc::setOption("ibm-shots","8192");
    auto buffer = acc->createBuffer("buffer", 2);

    auto compiler = xacc::getService<xacc::Compiler>("xacc-py");
    const std::string src = R"src(def f(buffer, t0):
       X(0)
       Ry(t0,1)
       CNOT(1,0)
       )src";

    auto ir = compiler->compile(src, acc);

    std::vector<double> p{-3.1415};
    auto f = ir->getKernel("f")->operator()(p);

    PauliOperator op;
    op.fromString("(5.9067,0) I + (-2.1433,0) X0 X1 + (-2.1433,0) Y0 Y1 + (.21829,0) Z0 + (-6.125,0) Z1");

    // std::cout << "OP: " << op.toString() << "\n";

    auto measureFunctions = op.toXACCIR()->getKernels();
    for (auto& m : measureFunctions) {
        m->insertInstruction(0,f);
    }

    PurificationDecorator decorator;
    decorator.setDecorated(acc);

    buffer->addExtraInfo("identity-coeff", ExtraInfo(5.9067));
    auto buffers = decorator.execute(buffer, measureFunctions);

    for (auto& b : buffers) b->print();
  }
}
int main(int argc, char **argv) {
  xacc::Initialize();
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
