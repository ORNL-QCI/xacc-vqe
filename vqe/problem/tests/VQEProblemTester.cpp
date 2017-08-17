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
#include <boost/math/constants/constants.hpp>

BOOST_AUTO_TEST_CASE(checkSimpleH2) {

    const std::string src = R"src(__qpu__ H2_sto-3g_singlet_H2_Molecule0_7414() {
	0.282342748105 1 1 0 1 0 0 1 0
	0.285028996894 2 1 0 1 0 0 2 0
	-0.941771530368 1 1 1 0
	0.285028996894 3 1 1 1 1 0 3 0
	0.297757778589 3 1 2 1 2 0 3 0
	0.285028996894 1 1 2 1 2 0 1 0
	-0.658507829418 2 1 2 0
	0.111545302727 2 1 1 1 3 0 0 0
	0.285028996894 1 1 3 1 3 0 1 0
	0.282342748105 0 1 1 1 1 0 0 0
	0.111545302727 3 1 2 1 0 0 1 0
	0.297757778589 2 1 3 1 3 0 2 0
	0.111545302727 2 1 0 1 2 0 0 0
	0.111545302727 3 1 0 1 2 0 1 0
	-0.941771530368 0 1 0 0
	0.111545302727 0 1 2 1 0 0 2 0
	0.111545302727 3 1 1 1 3 0 1 0
	-0.658507829418 3 1 3 0
	0.285028996894 3 1 0 1 0 0 3 0
	0.285028996894 0 1 2 1 2 0 0 0
	0.111545302727 0 1 3 1 1 0 2 0
	0.111545302727 1 1 2 1 0 0 3 0
	0.285028996894 2 1 1 1 1 0 2 0
	0.111545302727 0 1 1 1 3 0 2 0
	0.111545302727 1 1 3 1 1 0 3 0
	0.377702344336
	0.111545302727 2 1 3 1 1 0 0 0
	0.111545302727 1 1 0 1 2 0 3 0
	0.285028996894 0 1 3 1 3 0 0 0
})src";

	using namespace xacc::vqe;
	auto pi = boost::math::constants::pi<double>();

	xacc::setAccelerator("tnqvm");
	xacc::setOption("vqe-print-scaffold-source", "hello");
	xacc::Initialize(boost::unit_test::framework::master_test_suite().argc,
			boost::unit_test::framework::master_test_suite().argv);

	std::istringstream ss(src);

	VQEProblem<double> problem(ss);
	Eigen::VectorXd params(2);
	std::ofstream outFile("paramSweep.csv");

	params = -1.0*pi*Eigen::VectorXd::Ones(2) + (Eigen::VectorXd::Random(2)*0.5+Eigen::VectorXd::Ones(2)*0.5)*(pi-(-1*pi));
	std::cout << "Initial Params: " << params.transpose() << "\n";
//	params(0) = 3.14; //2.5;
//	params(1) = .21;
	std::cout << "HELLO: " << problem(params) << "\n";

	cppoptlib::NelderMeadSolver<VQEProblem<double>> solver;
	solver.setStopCriteria(VQEProblem<double>::getConvergenceCriteria());
	solver.minimize(problem, params);

	std::cout << "FINAL ENERGY: " << problem(params) << "\n";

//	for (auto angle1 = 2.5; angle1 <= pi; angle1 += .1) {
//		for (auto angle2 = 0.0; angle2 <= pi; angle2 += .1) {
//
//			params(0) = angle1;
//			params(1) = angle2;
//			auto energy = problem(params);
//			outFile << angle1 << ", " << angle2 << ", " << energy << "\n";
//			outFile.flush();
//		}
//	}

	outFile.close();

	xacc::Finalize();

}
