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
#include <gtest/gtest.h>
#include "DiagonalizeTask.hpp"
#include "ServiceRegistry.hpp"
#include "MPIProvider.hpp"

using namespace xacc::vqe;

TEST(ComputeEnergyVQETaskTester,checkSimple) {

	auto argc = xacc::getArgc();
	auto argv = xacc::getArgv();

	const std::string src = R"src(__qpu__ kernel() {
   0.7137758743754461
   -1.252477303982147 0 1 0 0
   0.337246551663004 0 1 1 1 1 0 0 0
   0.0906437679061661 0 1 1 1 3 0 2 0
   0.0906437679061661 0 1 2 1 0 0 2 0
   0.3317360224302783 0 1 2 1 2 0 0 0
   0.0906437679061661 0 1 3 1 1 0 2 0
   0.3317360224302783 0 1 3 1 3 0 0 0
   0.337246551663004 1 1 0 1 0 0 1 0
   0.0906437679061661 1 1 0 1 2 0 3 0
   -1.252477303982147 1 1 1 0
   0.0906437679061661 1 1 2 1 0 0 3 0
   0.3317360224302783 1 1 2 1 2 0 1 0
   0.0906437679061661 1 1 3 1 1 0 3 0
   0.3317360224302783 1 1 3 1 3 0 1 0
   0.3317360224302783 2 1 0 1 0 0 2 0
   0.0906437679061661 2 1 0 1 2 0 0 0
   0.3317360224302783 2 1 1 1 1 0 2 0
   0.0906437679061661 2 1 1 1 3 0 0 0
   -0.4759344611440753 2 1 2 0
   0.0906437679061661 2 1 3 1 1 0 0 0
   0.3486989747346679 2 1 3 1 3 0 2 0
   0.3317360224302783 3 1 0 1 0 0 3 0
   0.0906437679061661 3 1 0 1 2 0 1 0
   0.3317360224302783 3 1 1 1 1 0 3 0
   0.0906437679061661 3 1 1 1 3 0 1 0
   0.0906437679061661 3 1 2 1 0 0 1 0
   0.3486989747346679 3 1 2 1 2 0 3 0
   -0.4759344611440753 3 1 3 0
})src";

	std::shared_ptr<MPIProvider> provider;
	if (xacc::hasService<MPIProvider>("boost-mpi")) {
		provider = xacc::getService<MPIProvider>("boost-mpi");
	} else {
		provider = xacc::getService<MPIProvider>("no-mpi");
	}

	provider->initialize(argc,argv);
	auto world = provider->getCommunicator();
	xacc::setOption("n-qubits", "4");
	xacc::setOption("n-electrons", "2");
	xacc::setOption("vqe-task", "compute-energy");

	if (xacc::hasAccelerator("tnqvm")) {
		// Get the user-specified Accelerator,
		// or TNQVM if none specified
		auto accelerator = xacc::getAccelerator("tnqvm");

		DiagonalizeTask task;

		auto program = std::make_shared<VQEProgram>(accelerator, src, world);

		program->build();

		Eigen::VectorXd parameters(2);
		parameters << 0, 0;
		task.setVQEProgram(program);

		for (int i = 0; i < 10; i++) {
		   VQETaskResult result = task.execute(parameters);

		   EXPECT_NEAR(result.results[0].second, -1.13727042207, 1e-4);
		}
	}

}

int main(int argc, char** argv) {
   xacc::Initialize(argc,argv);
   ::testing::InitGoogleTest(&argc, argv);
   auto ret = RUN_ALL_TESTS();
   xacc::Finalize();
   return ret;
}
