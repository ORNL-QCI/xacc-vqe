#ifndef VQE_VQEPROBLEM_HPP_
#define VQE_VQEPROBLEM_HPP_

#include "problem.h"
#include "XACC.hpp"
#include <boost/math/constants/constants.hpp>
#include "GateQIR.hpp"
#include "IRGenerator.hpp"
#include "InstructionIterator.hpp"
#include "exprtk.hpp"
#include <boost/range/algorithm/count.hpp>

#include <boost/mpi.hpp>
#include "CountGatesOfTypeVisitor.hpp"

#include "solver/neldermeadsolver.h"
#include "solver/conjugatedgradientdescentsolver.h"
#include "solver/gradientdescentsolver.h"

using namespace xacc::quantum;

namespace xacc {

namespace vqe {

/**
 * The VQEProblem is a cppoptlib Problem subclass that
 * takes as input a file containing an XACC kernel
 * describing a molecular fermionic Hamiltonian. At
 * construction, the VQEProblem uses the FermionCompiler
 * to create a spin Hamiltonian representation of the
 * provided fermion Hamiltonian (via a Jordan Wigner, or other,
 * IR Transformation). This spin representation is composed
 * of XACC IR Functions, each representing a term in the
 * hamiltonian.
 *
 * VQEProblem provides a value() implementation that
 * computes the ground state energy of the input
 * molecule using N QPU invocations and computing
 * the corresponding expectation value of each. It computes
 * these expectation values in parallel and sums to produce
 * the ground state energy.
 *
 */
class VQEProblem : public cppoptlib::Problem<double> {

protected:

	using symbol_table_t = exprtk::symbol_table<double>;
	using expression_t = exprtk::expression<double>;
	using parser_t = exprtk::parser<double>;

	/**
	 * Reference to the QPU this VQE problem
	 * will run on.
	 */
	std::shared_ptr<Accelerator> qpu;

	/**
	 * The number of parameters in the
	 * state preparation circuit.
	 */
	int nParameters;

	/**
	 * The number of qubits used for this VQE problem.
	 */
	int nQubits;

	/**
	 * Reference to the state preparation circuit
	 * represented as XACC IR.
	 */
	std::shared_ptr<Function> statePrep;

	int rank = 0;

	int nRanks = 1;

	boost::mpi::communicator comm;

	void initialize(std::istream& moleculeKernel) {
		bool userProvidedKernels = false;

		// This class only takes kernels
		// represented as Fermion Kernels.
		xacc::setCompiler("fermion");

		// Create the Accelerator. This will be TNQVM
		// if --accelerator not passed to this executable.
		qpu = xacc::getAccelerator();

		// Count the number of kernels in the file
		std::string src((std::istreambuf_iterator<char>(moleculeKernel)),
		                 std::istreambuf_iterator<char>());

		auto nKernels = 0;
		size_t nPos = src.find("__qpu__", 0);
		while (nPos != std::string::npos) {
			nKernels++;
			nPos = src.find("__qpu__", nPos + 1);
		}

		std::vector<double> coeffs;
		// If nKernels > 1, we have non-fermioncompiler kernels
		// so lets check to see if they provided any coefficients
		if (nKernels > 1) { // && boost::contains(src, "coefficients")) {
			xacc::setCompiler("scaffold");
			userProvidedKernels = true;
		}

		// Create the Program
		Program program(qpu, src);
		program.addPreprocessor("fcidump-preprocessor");

		// Start compilation
		program.build();

		// Create a buffer of qubits
		nQubits = std::stoi(xacc::getOption("n-qubits"));

		// Get the Kernels that were created
		kernels = program.getRuntimeKernels();

		if (userProvidedKernels) {
			if (boost::contains(src, "pragma")
					&& boost::contains(src, "coefficient")) {
				std::vector<std::string> lines;
				boost::split(lines, src, boost::is_any_of("\n"));
				int counter = 0;
				for (int i = 0; i < lines.size(); ++i) {
					auto line = lines[i];
					if (boost::contains(line, "#pragma vqe-coefficient")) {
						std::vector<std::string> splitspaces;
						boost::split(splitspaces, line, boost::is_any_of(" "));
						boost::trim(splitspaces[2]);
						coeffs.push_back(std::stod(splitspaces[2]));
						InstructionParameter p(std::complex<double>(std::stod(splitspaces[2]), 0.0));
						InstructionParameter q((kernels[counter].getIRFunction()->nInstructions() == 0 ? 1 : 0));
						kernels[counter].getIRFunction()->addParameter(p);
						kernels[counter].getIRFunction()->addParameter(q);
						counter++;
					}
				}
			}
		}

		statePrep = createStatePreparationCircuit();

		// Set the number of VQE parameters
		nParameters = statePrep->nParameters();

		if (nParameters < 1) {
			XACCError("Error: State preparation circuit has 0 "
					"parameters. Try inputing a custom state prep kernel.");
		}

		if (xacc::optionExists("vqe-print-scaffold-source")) {
			printScaffoldSourceCode();
		}

		if (xacc::optionExists("vqe-print-stats")) {

			if (comm.rank() == 0) {
				XACCInfo("Number of Qubits = " + xacc::getOption("n-qubits"));
				XACCInfo("Number of Hamiltonian Terms = " + std::to_string(kernels.size()));
				XACCInfo("State Prep Type: " + statePrepType);
				XACCInfo("Number of Variational Parameters = " + std::to_string(nParameters));
				std::unordered_map<int, int> countKLocals;
				for (auto k : kernels) {
					// There is a measurement for each of the pauli's
					// in a given hamiltonian term, so just count them
					CountGatesOfTypeVisitor<Measure> visitor(k.getIRFunction());
					auto kLocalCount = visitor.countGates();
					auto search = countKLocals.find(kLocalCount);
					if (search != countKLocals.end()) {
						search->second++;
					} else {
						countKLocals[kLocalCount] = 1;
					}
				}

				std::map<int, int> ordered(countKLocals.begin(), countKLocals.end());
				for(auto& it: ordered) {
				     XACCInfo("N k-Local Terms (k,N) = ("
				    		 + std::to_string(it.first) + ", "
						 + std::to_string(it.second) + ")");
				}
			}

			exit(0);
		}
	}

public:

	/**
	 * Reference to the compiled XACC
	 * Kernels. These kernels each represent
	 * a term in the spin Hamiltonian (the
	 * Hamiltonian produced from a Jordan-Wigner
	 * or Bravi-Kitaev transformation). Each kernel
	 * amounts to a state preparation circuit
	 * followed by appropriately constructed qubit measurements.
	 */
	std::vector<Kernel<>> kernels;

	/**
	 * Referne
	 */
	int totalQpuCalls = 0;

	/**
	 * Reference to the current iteration's computed energy.
	 */
	double currentEnergy;

	VQEProblem(std::istream& moleculeKernel, boost::mpi::communicator& c) :
			comm(c), nParameters(0), currentEnergy(0.0) {
		initialize(moleculeKernel);
	}

	/**
	 * The constructor, executes the compilation from
	 * a fermionic hamiltonian represented as an XACC kernel
	 * function and provided as a file stream.
	 * Compilation maps this hamiltonian to a spin-based
	 * hamiltonian through a Jordan-Wigner (or other) IR Transformation. This
	 * constructor also prepares a user-specified state preparation circuit.
	 *
	 * @param moleculeKernel
	 */
	VQEProblem(std::istream& moleculeKernel) :
			nParameters(0), currentEnergy(0.0), comm(boost::mpi::communicator()) {
		initialize(moleculeKernel);
	}

	virtual const Eigen::VectorXd minimize() {
		cppoptlib::NelderMeadSolver<VQEProblem> solver;
		solver.setStopCriteria(VQEProblem::getConvergenceCriteria());
		auto params = initializeParameters();
		solver.minimize(*this, params);
		return params;
	}

	/**
	 * Return an initial random vector of
	 * parameters.
	 *
	 * @return params Initial parameters
	 */
	Eigen::VectorXd initializeParameters() {
		std::srand(time(0));
		auto pi = boost::math::constants::pi<double>();
		Eigen::VectorXd rand;

		std::vector<double> data;

		// Random parameters between -pi and pi
		if (comm.rank() == 0) {
			rand = -1.0 * pi * Eigen::VectorXd::Ones(nParameters)
					+ (Eigen::VectorXd::Random(nParameters) * 0.5
							+ Eigen::VectorXd::Ones(nParameters) * 0.5)
							* (pi - (-1 * pi));
			data.resize(rand.size());
			Eigen::VectorXd::Map(&data[0], rand.size()) = rand;
		}

		boost::mpi::broadcast(comm, data, 0);

		if (comm.rank() != 0) {
			rand = Eigen::Map<Eigen::VectorXd>(data.data(), data.size());
		}

		return rand;
	}

	/**
	 * Compute the energy at the provided vector
	 * of initial state parameters.
	 *
	 * @param x The parameters
	 * @return energy The computed energy.
	 */
	double value(const Eigen::VectorXd& x) {
		// Local Declarations
		double sum = 0.0, localExpectationValue = 0.0;

		// Evaluate our variable parameterized State Prep circuite
		// to produce a state prep circuit with actual rotations
		auto evaluatedStatePrep = evaluateStatePreparationCircuit(x);

		// Execute the kernels on the appropriate QPU
		// in parallel using OpenMP threads per
		// every MPI rank.
		int rank = comm.rank();
		int nRanks = comm.size();
		int myStart = (rank) * kernels.size() / nRanks;
		int myEnd = (rank + 1) * kernels.size() / nRanks;
#pragma omp parallel for reduction (+:sum)
		for (int i = myStart; i < myEnd; i++) {

			// Get the ith Kernel
			auto kernel = kernels[i];

			// Insert the state preparation circuit IR
			// at location 0 in this Kernels IR instructions.
			kernel.getIRFunction()->insertInstruction(0, evaluatedStatePrep);

			// Create a temporary buffer of qubits
			auto buff = qpu->createBuffer("qreg", nQubits);

			if (kernel.getIRFunction()->nInstructions() > 1) {
				// Execute the kernel!
				kernel(buff);
			}

			// Get Expectation value. The second parameter of
			// the Kernel's IR function stores whether or not
			// this kernel is the Identity.
			if (boost::get<int>(kernel.getIRFunction()->getParameter(1))) {
				localExpectationValue = 1.0;
			} else {
				localExpectationValue = buff->getExpectationValueZ();
			}

			// Sum up the expectation values, the Hamiltonian
			// terms coefficient is stored in the first
			// parameter of the Kernels IR Function representation
			sum += std::real(
					boost::get<std::complex<double>>(
							kernel.getIRFunction()->getParameter(0)))
					* localExpectationValue;

			// The next iteration will have a different
			// state prep circuit, so toss the current one.
			kernel.getIRFunction()->removeInstruction(0);
		}

		// Set the energy.
		currentEnergy = sum;

		double result = 0.0;
		boost::mpi::all_reduce(comm, currentEnergy, result, std::plus<double>());
		currentEnergy = result;

		totalQpuCalls += kernels.size();

		std::stringstream ss;
		ss << x.transpose();

		if (rank == 0) XACCInfo("Computed VQE Energy = " + std::to_string(currentEnergy) + " at (" + ss.str() + ")");
		return currentEnergy;
	}

	/**
	 * This internal structure describes the criteria for
	 * convergence in this VQE execution.
	 */
	class VQECriteria : public cppoptlib::Criteria<double> {
	public:
	    static VQECriteria defaults() {
	    	VQECriteria d;
	        d.iterations = 1000;
	        d.xDelta = 0;
	        d.fDelta = 1e-6;
	        d.gradNorm = 1e-4;
	        d.condition = 0;
	        return d;
	    }
	};

	/**
	 * This static method returns the user-specified
	 * convergence criteria.
	 *
	 * @return criteria The convergence criteria.
	 */
	static VQECriteria getConvergenceCriteria() {
		auto criteria = VQECriteria::defaults();

		if (xacc::optionExists("vqe-energy-delta")) {
			criteria.fDelta = std::stod(xacc::getOption("vqe-energy-delta"));
		}

		if (xacc::optionExists("vqe-iterations")) {
			criteria.iterations = std::stoi(xacc::getOption("vqe-iterations"));
		}

		return criteria;

	}

	const int getNParameters() {
		return nParameters;
	}

private:

	std::string statePrepType = "uccsd";

	/**
	 * This private utility method is invoked once to
	 * produce the state preparation circuit, represented
	 * as XACC IR.
	 *
	 * @return statePrep The Function representation of
	 * 			the State Preparation circuite
	 */
	std::shared_ptr<Function> createStatePreparationCircuit() {

		if (xacc::optionExists("vqe-state-prep-kernel")) {
			auto filename = xacc::getOption("vqe-state-prep-kernel");
			std::ifstream filess(filename);

			xacc::setCompiler("scaffold");
			if (xacc::optionExists("vqe-state-prep-kernel-compiler")) {
				xacc::setCompiler(xacc::getOption("vqe-state-prep-kernel-compiler"));
			}

			statePrepType = "custom";

			Program p(qpu, filess);
			p.build();

			auto kernel = p.getRuntimeKernels()[0];

			return kernel.getIRFunction();
		} else {
			if (xacc::optionExists("state-preparation")) {
				statePrepType = xacc::getOption("state-preparation");
			}

			auto statePrepGenerator = ServiceRegistry::instance()->getService<
					IRGenerator>(statePrepType);
			return statePrepGenerator->generate(std::make_shared<AcceleratorBuffer>("",nQubits));
		}
	}

	/**
	 * This private utility method is used to return a
	 * new IR Function representing the state preparation circuit
	 * evaluated at the given parameter vector.
	 *
	 * @param x The parameters
	 * @return evaluatedFunction Evaluated State Preparation circuit
	 */
	std::shared_ptr<Function> evaluateStatePreparationCircuit(const Eigen::VectorXd& x) {

		auto pi = boost::math::constants::pi<double>();
		std::vector<std::string> variableNames;
		for (int i = 0; i < x.rows(); i++) {
			variableNames.push_back(
					boost::get<std::string>(statePrep->getParameter(i)));
		}

		auto evaluatedStatePrep = std::make_shared<GateFunction>("stateprep");
		for (auto inst : statePrep->getInstructions()) {
			if (!inst->isComposite() && inst->isParameterized()
					&& inst->getParameter(0).which() == 3) {
				int idx = -1;
				auto expression = boost::get<std::string>(inst->getParameter(0));
				for (int i = 0; i < nParameters; i++) {
					if (boost::contains(expression, variableNames[i])) {
						idx = i;
					}
				}
				std::string varName = variableNames[idx];
				double val;
				symbol_table_t symbol_table;
				symbol_table.add_variable(varName, val);
				symbol_table.add_constants();
				expression_t expr;
				expr.register_symbol_table(symbol_table);
				parser_t parser;
				parser.compile(expression, expr);
				val = x(idx);
				auto res = expr.value();
				if (res < 0.0) {
					res = 4 * pi + res;
				}
				InstructionParameter p(res);
				auto updatedInst = GateInstructionRegistry::instance()->create(
						inst->getName(), inst->bits());
				updatedInst->setParameter(0, p);
				evaluatedStatePrep->addInstruction(updatedInst);
			} else {
				evaluatedStatePrep->addInstruction(inst);
			}
		}

		return evaluatedStatePrep;

	}

	/**
	 * This method prints the compiled Hamiltonian kernels
	 * to Scaffold source code.
	 */
	void printScaffoldSourceCode() {
		int counter = 0;
		for (auto k : kernels) {
			auto scaffold = xacc::getCompiler("scaffold");
			auto f = k.getIRFunction();
			f->insertInstruction(0, statePrep);
			auto srcStr = scaffold->translate("qreg", f);
			f->removeInstruction(0);

			boost::filesystem::path dir("scaffold_source");
			if (!boost::filesystem::exists(dir)) {
				if (!boost::filesystem::create_directory(dir)) {
					XACCError("Could not create scaffold_source directory.");
				}
			}
			std::ofstream out(
					"scaffold_source/"
							+ xacc::getOption("vqe-print-scaffold-source")
							+ "_kernel_" + std::to_string(counter) + ".hpp");
			out << srcStr;
			out.flush();
			out.close();
			counter++;
		}

		if (xacc::optionExists("vqe-exit-after-scaffold")) {
			xacc::Finalize();
			exit(0);
		}
	}
};
}

}

#endif
//
//auto scaffold = xacc::getCompiler("scaffold");
//	auto f = kernel.getIRFunction();
//	auto srcStr = scaffold->translate("qreg", f);
//
//	boost::filesystem::path dir("temp");
//	if (!boost::filesystem::exists(dir)) {
//		if (!boost::filesystem::create_directory(dir)) {
//			XACCError("Could not create scaffold_source directory.");
//		}
//	}
//	std::ofstream out("temp/kernel_" + std::to_string(i) + "_at_"
//					+ std::to_string(boost::get<double>(parameters[0])) + "_" + std::to_string(boost::get<double>(parameters[1]))
//					+ ".hpp");
//	out << srcStr;
//	out.flush();
//	out.close();

/*
			auto statePrep = std::dynamic_pointer_cast<Function>(k.getIRFunction()->getInstruction(0));
			InstructionIterator it2(statePrep);
			it2.next();
			int counter = 0, index = 0;
			while (it2.hasNext()) {
				auto next = it2.next();
				if (next->getName() == "Rz"
						&& std::fabs(boost::get<double>(next->getParameter(0)))
								< 1e-12) {
//					std::cout << "Found an Rz with 0 angle\n";

//					next->disable();
					// If I'm here, I have counter instructions in front of me
					// So disable the counter in front and counter behind
					while(counter != 0) {
//					for (int i = index-counter; i <= 2*(counter+1); i++) {
//						std::cout << counter << ", " << index << " | Disabling " << statePrep->getInstruction(index-counter)->toString("qreg") << " and " << statePrep->getInstruction(index+counter)->toString("qreg") << "\n";
//						statePrep->getInstruction(index-counter)->disable();
//						statePrep->getInstruction(index+counter)->disable();
						counter--;
					}

					// reset the counter
//					counter = 0;
				}

				if (next->getName() != "X" && next->isEnabled()) {
					counter++;
				}

				index++;
			}*/
