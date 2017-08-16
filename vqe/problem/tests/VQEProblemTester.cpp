/*
 * VQEProblemTester.hpp
 *
 *  Created on: Aug 14, 2017
 *      Author: aqw
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE VQEProblemTester

#include <boost/test/included/unit_test.hpp>
#include "VQEProblem.hpp"
#include "solver/neldermeadsolver.h"

BOOST_AUTO_TEST_CASE(checkSimpleH2) {

    const std::string src = R"src(__qpu__ H2_sto-3g_singlet_H2_Molecule0_1645885() {
	0.37524014544 1 1 0 1 0 0 1 0
	0.370097420049 2 1 0 1 0 0 2 0
	-1.54693295173 1 1 1 0
	0.370097420049 3 1 1 1 1 0 3 0
	0.391450513269 3 1 2 1 2 0 3 0
	0.370097420049 1 1 2 1 2 0 1 0
	0.0268523970537 2 1 2 0
	0.0805940580404 2 1 1 1 3 0 0 0
	0.370097420049 1 1 3 1 3 0 1 0
	0.37524014544 0 1 1 1 1 0 0 0
	0.0805940580404 3 1 2 1 0 0 1 0
	0.391450513269 2 1 3 1 3 0 2 0
	0.0805940580404 2 1 0 1 2 0 0 0
	0.0805940580404 3 1 0 1 2 0 1 0
	-1.54693295173 0 1 0 0
	0.0805940580404 0 1 2 1 0 0 2 0
	0.0805940580404 3 1 1 1 3 0 1 0
	0.0268523970537 3 1 3 0
	0.370097420049 3 1 0 1 0 0 3 0
	0.370097420049 0 1 2 1 2 0 0 0
	0.0805940580404 0 1 3 1 1 0 2 0
	0.0805940580404 1 1 2 1 0 0 3 0
	0.370097420049 2 1 1 1 1 0 2 0
	0.0805940580404 0 1 1 1 3 0 2 0
	0.0805940580404 1 1 3 1 1 0 3 0
	1.70138568668
	0.0805940580404 2 1 3 1 1 0 0 0
	0.0805940580404 1 1 0 1 2 0 3 0
	0.370097420049 0 1 3 1 3 0 0 0
)src";

	using namespace xacc::vqe;

	xacc::Initialize(boost::unit_test::framework::master_test_suite().argc,
			boost::unit_test::framework::master_test_suite().argv);

	std::istringstream ss(src);

	VQEProblem<double> problem(ss);
	Eigen::VectorXd params(2);
	params << 0, 0.05677;

	std::cout << "ENERGY: " << problem(params) << "\n";

	xacc::Finalize();

}
