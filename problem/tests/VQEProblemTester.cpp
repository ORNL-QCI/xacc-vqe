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
#include <boost/mpi.hpp>

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

    const std::string src02 = R"src02(__qpu__ H2_sto-3g_singlet_H2_Molecule0_2() {
        0.370213650434 1 1 0 1 0 0 1 0
        0.36458386954 2 1 0 1 0 0 2 0
        -1.49816537423 1 1 1 0
        0.36458386954 3 1 1 1 1 0 3 0
        0.384950478934 3 1 2 1 2 0 3 0
        0.36458386954 1 1 2 1 2 0 1 0
        -0.0843846792728 2 1 2 0
        0.0818185944671 2 1 1 1 3 0 0 0
        0.36458386954 1 1 3 1 3 0 1 0
        0.370213650434 0 1 1 1 1 0 0 0
        0.0818185944671 3 1 2 1 0 0 1 0
        0.384950478934 2 1 3 1 3 0 2 0
        0.0818185944671 2 1 0 1 2 0 0 0
        0.0818185944671 3 1 0 1 2 0 1 0
        -1.49816537423 0 1 0 0
        0.0818185944671 0 1 2 1 0 0 2 0
        0.0818185944671 3 1 1 1 3 0 1 0
        -0.0843846792728 3 1 3 0
        0.36458386954 3 1 0 1 0 0 3 0
        0.36458386954 0 1 2 1 2 0 0 0
        0.0818185944671 0 1 3 1 1 0 2 0
        0.0818185944671 1 1 2 1 0 0 3 0
        0.36458386954 2 1 1 1 1 0 2 0
        0.0818185944671 0 1 1 1 3 0 2 0
        0.0818185944671 1 1 3 1 1 0 3 0
        1.40014259045
        0.0818185944671 2 1 3 1 1 0 0 0
        0.0818185944671 1 1 0 1 2 0 3 0
        0.36458386954 0 1 3 1 3 0 0 0
})src02";

	using namespace xacc::vqe;
	auto pi = boost::math::constants::pi<double>();

	xacc::setAccelerator("tnqvm");
	xacc::setOption("vqe-print-scaffold-source", "hello");
	xacc::setOption("n-electrons", "2");

	xacc::Initialize(boost::unit_test::framework::master_test_suite().argc,
			boost::unit_test::framework::master_test_suite().argv);

	boost::mpi::environment env;
	boost::mpi::communicator comm;

	std::istringstream ss(src);

	VQEProblem problem(ss, comm);
	Eigen::VectorXd params(2);
	std::ofstream outFile("paramSweep.csv");

	params(0) = 1.57; //2.5;
	params(1) = -.75;
//	cppoptlib::NelderMeadSolver<VQEProblem> solver;
//	solver.setStopCriteria(VQEProblem::getConvergenceCriteria());
//
//	solver.minimize(problem, params);

	std::cout << "FINAL ENERGY: " << problem(params) << "\n";
/*
	Eigen::VectorXd range = Eigen::VectorXd::LinSpaced(100, -pi, pi);
	std::vector<double> rangeVec(range.data(),
			range.data() + range.size());
	for (auto theta0 : rangeVec) {
		for (auto theta1 : rangeVec) {
			params(0) = theta0;
			params(1) = theta1;
			outFile << theta0 << ", " << theta1 << ", " << problem(params)
					<< "\n";
			outFile.flush();
		}
	}
*/
	outFile.close();

	xacc::Finalize();

}
//
