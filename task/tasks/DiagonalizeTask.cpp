#include "FermionToSpinTransformation.hpp"
#include <iostream>
#include <unordered_map>

#include "DiagonalizeTask.hpp"
#include <Eigen/Sparse>

namespace xacc {
namespace vqe {


VQETaskResult DiagonalizeTask::execute(
		Eigen::VectorXd parameters) {
	int nQubits = std::stoi(xacc::getOption("n-qubits"));

	auto hamiltonianInstruction = program->getPauliOperator();

	std::shared_ptr<DiagonalizeBackend> backend;
	if (xacc::optionExists("diagonalize-backend")) {
		auto str = xacc::getOption("diagonalize-backend");
		backend = ServiceRegistry::instance()->getService<
				DiagonalizeBackend>(str);
	} else {
		backend = ServiceRegistry::instance()->getService<
				DiagonalizeBackend>("diagonalize-eigen");
	}

	auto energy = backend->diagonalize(program);

	VQETaskResult result;
	result.angles = parameters;
	result.energy = energy;

	result.results.push_back({parameters, energy});
	return result;
}

double EigenDiagonalizeBackend::diagonalize(
		std::shared_ptr<VQEProgram> prog) {
	std::complex<double> gsReal;
	auto hamiltonian = prog->getPauliOperator();
//	auto nTerms = inst.nTerms();
	auto nQubits = prog->getNQubits();

	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	auto getBitStrForIdx = [&](std::uint64_t i) {
		std::stringstream s;
		for (int k = nQubits - 1; k >= 0; k--) s << ((i >> k) & 1);
		return s.str();
	};

	Eigen::VectorXd eigenvalues;
	if (xacc::optionExists("n-electrons")) {
		int nElectrons = std::stoi(xacc::getOption("n-electrons"));

		// Generate all n-qubit bitstrings with n-electron
		// bits set
		std::string initBitString = "";
		for (int i = 0; i < nQubits - nElectrons; i++)
			initBitString += "0";
		for (int i = 0; i < nElectrons; i++)
			initBitString += "1";

//		std::cout << "Initial bit string = " << initBitString << "\n";

		std::vector<std::string> bitStrings;
		do {
//			std::cout << initBitString << "\n";
			bitStrings.push_back(initBitString);
		} while (std::next_permutation(initBitString.begin(),
				initBitString.end()));

		std::map<std::string, std::uint64_t> bitsToIdx;
		std::map<std::uint64_t, std::string> idxToBits;

		std::uint64_t counter = 0;
		for (auto& bs : bitStrings) {
			bitsToIdx.insert( { bs, counter });
			idxToBits.insert( { counter, bs });
			counter++;
		}
		int nBitStrings = counter;

		xacc::info("Considering Hamiltonian subspace spanned by "
				+ std::to_string(nBitStrings) + " eigenstates with "
				+ std::to_string(nElectrons) + " occupations");

		Eigen::MatrixXcd mat(nBitStrings, nBitStrings);
		mat.setZero();
		for (std::uint64_t i = 0; i < nBitStrings; i++) {
			auto bitStr = idxToBits[i];
			auto braResults = hamiltonian.computeActionOnKet(bitStr);

			for (auto& r : braResults) {
				if (bitStr == r.first) {
					mat(i, i) += r.second;
				}
			}

			for (auto& result : braResults) {
				std::uint64_t k = bitsToIdx[result.first];
				if (i != k)
					mat(i, k) += result.second;
			}

		}
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(mat);

		eigenvalues = es.eigenvalues();
	} else {
		Eigen::MatrixXcd A(dim, dim);
		for (std::uint64_t myRow = 0; myRow < dim; myRow++) {
			auto rowBitStr = getBitStrForIdx(myRow);
			auto results = hamiltonian.computeActionOnBra(rowBitStr);
			for (auto& result : results) {
				std::uint64_t k = std::stol(result.first, nullptr, 2);
				A(myRow, k) += result.second;
			}
		}
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);

		eigenvalues = es.eigenvalues();
	}
	std::cout << "Spectrum: " << eigenvalues.transpose() << "\n";
	gsReal = eigenvalues(0);
	std::stringstream ss;
	ss << std::setprecision(12) << gsReal;
	xacc::info("Ground State Energy of Hamiltonian = " + ss.str());
	return std::real(gsReal);
}

}
}
