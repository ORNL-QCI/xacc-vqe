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
#include "Eigen/Dense"

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


  PauliOperator expected, actual, op;
  op.fromXACCIR(newIR);
  expected.fromString("(-0.335683,0) I + (-0.780643,0) Z0 + (0.181625,0) X0");

  actual.fromXACCIR(newIR);

  EXPECT_TRUE(actual == expected);

   auto data = op.toDenseMatrix(6).data();
   Eigen::MatrixXcd A = Eigen::Map<Eigen::MatrixXcd>(data, 64,64);
//    auto A = op.toDenseMatrix(6);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
    auto reducedEnergy = es.eigenvalues()[0];

    std::cout << "SHOULD BE: " << std::setprecision(12) << reducedEnergy << "\n";
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

    std::vector<double> p{-3.1415,1.57};
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

TEST(QubitTaperingTester, checkSixQubit) {

  if (xacc::hasAccelerator("local-ibm")) {
    auto acc = xacc::getAccelerator("local-ibm");
    auto b = acc->createBuffer("q",4);

    auto compiler = xacc::getService<xacc::Compiler>("xacc-py");
    const std::string src2 = R"src(def f(buffer, t0, t1):
       Ry(t0,0)
       Rz(t1,0)
       )src";

    auto ir = compiler->compile(src2, acc);

    std::vector<double> p{-3.1415,1.57};
    auto f = ir->getKernel("f")->operator()(p);

    PauliOperator op;
    op.fromString("(-0.00421023,0) Y3 Y4 Z5 + (0.0816923,0) Z3 Z5 + (-0.0140516,0) Y3 Y5 + (-0.00895403,0) Z3 Y4 Y5 + (-0.00369882,0) Z2 X4 X5 + (-0.00369882,0) Z2 Y4 Y5 + (-0.00944179,0) Z2 Y3 Z4 Y5 + (0.101934,0) Z2 Z3 + (-0.00369882,0) X1 X2 Z5 + (-0.00369882,0) Y1 Y2 Z5 + (0.0165034,0) X1 X2 X4 X5 + (0.0165034,0) Y1 Y2 X4 X5 + (0.0165034,0) Y1 Y2 Y4 Y5 + (-0.00895403,0) Z3 X4 X5 + (0.0138874,0) Z1 X4 X5 + (-0.00277645,0) X1 X2 X3 Z4 X5 + (-0.00421023,0) X3 X4 Z5 + (-0.00277645,0) X1 X2 Y3 Z4 Y5 + (-0.00277645,0) Y1 Y2 X3 Z4 X5 + (-0.00677708,0) X1 X2 X3 X4 + (-0.00677708,0) X1 X2 Y3 Y4 + (-0.00677708,0) Y1 Y2 X3 X4 + (-0.00677708,0) Y1 Y2 Y3 Y4 + (-0.0218996,0) Y1 Y2 Z3 + (-0.00421023,0) X0 X1 Z2 + (0.0129455,0) X0 X1 Y3 Z4 Y5 + (0.0112333,0) Y0 Y1 X3 X4 + (0.0821053,0) Z4 Z5 + (0.0254908,0) Y0 Y1 Z3 + (0.0709831,0) Z0 Z4 + (-0.0140516,0) X3 X5 + (-0.00677708,0) Y0 Y1 X4 X5 + (-0.00677708,0) X0 X1 Y4 Y5 + (0.0229208,0) Z0 Y3 Z4 Y5 + (-0.00895403,0) Z0 Y1 Y2 + (0.11747,0) Z2 Z5 + (0.0129455,0) Y0 Z1 Y2 Y3 Y4 + (0.0138874,0) X1 X2 Z4 + (-0.00677708,0) X0 X1 X4 X5 + (0.0254908,0) Z0 Y3 Y4 + (0.158901,0) Z0 Z3 + (-0.00421023,0) Y0 Y1 Z2 + (-0.00944179,0) Y0 Z1 Y2 Z5 + (0.110544,0) Z1 Z4 + (-0.00883391,0) Z1 Y3 Y4 + (0.0254908,0) Z0 X3 X4 + (-0.115794,0) Z0 + (-0.0140516,0) X0 X2 + (-0.00895403,0) Z0 X1 X2 + (0.0821053,0) Z1 Z2 + (0.0597498,0) Z0 Z1 + (0.0138874,0) Z1 Y4 Y5 + (-0.568532,0) Z5 + (-0.0410422,0) X4 X5 + (0.0536886,0) X0 Z1 X2 + (0.0816923,0) Z0 Z2 + (0.0399687,0) X3 X4 + (0.0536886,0) Y0 Z1 Y2 + (0.0597498,0) Z3 Z4 + (-158.749,0) + (0.0399687,0) X0 X1 + (0.0399687,0) Y0 Y1 + (0.0399687,0) Y3 Y4 + (-0.0410422,0) Y4 Y5 + (-0.00143379,0) Z2 Y3 Y4 + (-0.382003,0) Z1 + (-0.0410422,0) Y1 Y2 + (0.0229208,0) Z0 X3 Z4 X5 + (-0.115794,0) Z3 + (0.0112333,0) Y0 Y1 Y3 Y4 + (0.0536886,0) Y3 Z4 Y5 + (0.0165034,0) X1 X2 Y4 Y5 + (-0.00277645,0) Y1 Y2 Y3 Z4 Y5 + (-0.0140516,0) Y0 Y2 + (-0.568532,0) Z2 + (0.0112333,0) X0 X1 X3 X4 + (0.0986087,0) Z2 Z4 + (0.0536886,0) X3 Z4 X5 + (0.0129455,0) X0 X1 X3 Z4 X5 + (-0.382003,0) Z4 + (0.0112333,0) X0 X1 Y3 Y4 + (0.0129455,0) Y0 Y1 Y3 Z4 Y5 + (0.0202421,0) X0 Z1 X2 Y3 Z4 Y5 + (0.0129455,0) Y0 Y1 X3 Z4 X5 + (0.0229208,0) Y0 Z1 Y2 Z3 + (0.0709831,0) Z1 Z3 + (-0.0072745,0) Y0 Z1 Y2 Z4 + (0.0229208,0) X0 Z1 X2 Z3 + (-0.0218996,0) Z0 Y4 Y5 + (-0.00944179,0) Z2 X3 Z4 X5 + (0.0129455,0) Y0 Z1 Y2 X3 X4 + (0.0129455,0) X0 Z1 X2 Y3 Y4 + (0.0129455,0) X0 Z1 X2 X3 X4 + (0.0138874,0) Y1 Y2 Z4 + (-0.0072745,0) X0 Z1 X2 Z4 + (0.101934,0) Z0 Z5 + (0.0202421,0) Y0 Z1 Y2 X3 Z4 X5 + (0.0254908,0) X0 X1 Z3 + (0.0202421,0) X0 Z1 X2 X3 Z4 X5 + (-0.00143379,0) Y0 Y1 Z5 + (-0.0218996,0) X1 X2 Z3 + (-0.0072745,0) Z1 Y3 Z4 Y5 + (-0.00277645,0) Y0 Z1 Y2 Y4 Y5 + (-0.00883391,0) Y0 Y1 Z4 + (-0.00883391,0) X0 X1 Z4 + (-0.00143379,0) Z2 X3 X4 + (-0.0218996,0) Z0 X4 X5 + (-0.00677708,0) Y0 Y1 Y4 Y5 + (-0.00277645,0) Y0 Z1 Y2 X4 X5 + (-0.00277645,0) X0 Z1 X2 Y4 Y5 + (-0.00277645,0) X0 Z1 X2 X4 X5 + (0.0986087,0) Z1 Z5 + (-0.00143379,0) X0 X1 Z5 + (-0.0410422,0) X1 X2 + (0.0202421,0) Y0 Z1 Y2 Y3 Z4 Y5 + (-0.00944179,0) X0 Z1 X2 Z5 + (-0.00883391,0) Z1 X3 X4 + (-0.0072745,0) Z1 X3 Z4 X5");
    auto tir = op.toXACCIR();

    QubitTapering tapering;
    auto newIR = tapering.transform(tir);

    PauliOperator expected, actual;
    expected.fromString("(0.0072745,-0) Z0 Z2 X3 + (0.0229208,0) Z0 X1 Z2 Z3 + (-0.101934,-0) Z0 Z1 Z3 + (0.0138874,0) X0 X1 Z2 + (-0.00369882,0) Y0 Y1 Z3 + (-0.628282,0) Z3 + (-0.0709831,-0) Z0 Z2 Z3 + (-0.0499962,0) Y2 Y3 + (0.00277645,-0) Z0 X1 Y2 Y3 + (0.00677708,-0) X0 X2 X3 + (-0.0499962,0) X0 X1 + (0.00277645,-0) X0 X1 Z2 X3 + (0.0218996,-0) Z0 Z1 Y2 Y3 + (0.0218996,-0) X0 X1 Z2 Z3 + (0.0112333,0) X0 X2 + (-0.00369882,0) X0 X1 Z3 + (0.00677708,0) Y0 Y1 X2 Z3 + (0.00277645,0) X1 Y2 Y3 + (-0.101934,-0) Z1 Z2 Z3 + (-0.0357585,0) X2 Z3 + (0.158901,0) Z0 Z1 Z2 Z3 + (0.00277645,-0) Y0 Y1 Z2 X3 + (0.0254908,0) X0 Z2 Z3 + (-0.039637,0) Z2 X3 + (-0.039637,0) X1 + (0.0254908,-0) Z0 Z1 X2 Z3 + (0.0165034,0) Y0 Y1 Y2 Y3 + (-0.0499962,0) X2 X3 + (0.110544,0) Z0 Z2 + (0.0202421,0) X1 X3 + (0.00883391,-0) Z0 X2 + (0.0072745,0) Z0 X3 + (0.0986087,0) Z1 Z2 + (0.0254908,-0) X0 Z1 Z2 Z3 + (0.0986087,0) Z0 Z3 +(0.0072745,0) X1 Z2 + (-0.00369882,0) Z1 X2 X3 + (-0.00369882,0) Z1 Y2 Y3 + (0.0129455,-0) X0 Z1 Z2 X3 + (0.0202421,0) Z0 X1 Z2 X3 + (0.00277645,-0) Z0 X1 X2 X3 + (-0.0499962,0) Y0 Y1 + (0.0202421,-0) X1 Z2 X3 +(0.00143379,0) X0 Z1 Z3 + (-0.0357585,0) X2 + (0.00677708,-0) X0 Y2 Y3 + (-0.628282,0) Z1 + (0.0138874,0)Z0 X2 X3 + (-0.0357585,0) X0 Z1 + (-0.039637,0) Z0 X1 + (0.00677708,0) X0 X1 X2 Z3 + (-0.463695,0) Z2 + (0.0165034,0) X0 X1 X2 X3 + (0.00143379,-0) Z1 X2 + (0.0138874,0) Z0 Y2 Y3 + (0.0112333,-0) X0 X2 Z3 + (0.0138874,0) Y0 Y1 Z2 + (0.00883391,-0) X0 Z2 + (0.00677708,0) X0 Z1 Y2 Y3 + (0.0165034,0) Y0 Y1 X2 X3 + (0.00944179,0) Z1 X3 + (0.0129455,-0) Z0 X1 X2 Z3 + (0.0129455,-0) X0 X3 + (0.0229208,0) Z0 Z1 Z2 X3 + (0.00143379,-0) X0 Z3 + (0.0254908,0) Z0 Z1 X2 + (0.00944179,-0) Z0 X1 Z3 + (-158.749,0) + (0.0218996,-0) Z0 Z1 X2 X3 + (0.0112333,0) X0 Z1 X2 Z3 + (0.197899,0) Z0 Z1 + (0.00944179,-0) Z1 Z2 X3 + (0.00883391,0) Z0 X2 Z3+ (-0.039637,0) X3 + (0.0218996,-0) Y0 Y1 Z2 Z3 + (0.0129455,-0) X1 X2 + (0.00277645,0) X1 X2 X3 + (0.0202421,-0) Z0 X1 X3 + (0.197899,0) Z2 Z3 + (0.0229208,-0) X1 Z2 Z3 + (-0.0709831,-0) Z0 Z1 Z2 + (0.00277645,0) X0 X1 X3 + (0.0112333,-0) X0 Z1 X2 + (0.00944179,0) X1 Z3 + (0.00277645,0) Y0 Y1 X3 + (0.0129455,0) X0 Z1 X3 + (0.0129455,0) X1 X2 Z3 + (0.11747,0) Z1 Z3 + (0.0072745,-0) Z0 X1 Z2 + (-0.0357585,0) X0 + (0.00677708,0) X0 Z1 X2 X3 + (0.00143379,0) Z1 X2 Z3 + (0.0129455,0) Z0 X1 X2 + (0.00883391,0) X0 Z1 Z2 + (0.0165034,0) X0 X1 Y2 Y3 + (0.00677708,-0) Y0 Y1 X2 + (0.0129455,0) X0 Z2 X3 + (0.0229208,-0) Z0 Z1 X3 + (0.00677708,-0) X0 X1 X2 + (-0.463695,0) Z0");

    actual.fromXACCIR(newIR);

    EXPECT_TRUE(4 == actual.nQubits());

    auto data = op.toDenseMatrix(6).data();
    Eigen::MatrixXcd A = Eigen::Map<Eigen::MatrixXcd>(data, 64,64);

    // auto A = op.toDenseMatrix(6);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
    auto reducedEnergy = es.eigenvalues()[0];

    auto datab = actual.toDenseMatrix(4).data();
    Eigen::MatrixXcd B = Eigen::Map<Eigen::MatrixXcd>(datab, 16,16);
    // auto B = actual.toDenseMatrix(4);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es2(B);
    auto reducedEnergy2 = es2.eigenvalues()[0];
    std::cout << "SHOULD BE: " << std::setprecision(12) << reducedEnergy << "\n";

    // EXPECT_NEAR(reducedEnergy, reducedEnergy2, 1e-5);

  }
}
int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
