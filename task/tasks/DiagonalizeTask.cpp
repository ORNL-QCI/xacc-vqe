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
		backend = xacc::getService<
				DiagonalizeBackend>(str);
	} else {
		backend = xacc::getService<
				DiagonalizeBackend>("diagonalize-eigen");
	}

	auto energy = backend->diagonalize(program);

	VQETaskResult result;
	result.angles = parameters;
	result.energy = energy;
	return result;
}

double EigenDiagonalizeBackend::diagonalize(
		std::shared_ptr<VQEProgram> prog) {
	std::complex<double> gsReal;
	auto hamiltonian = prog->getPauliOperator();
	auto nQubits = prog->getNQubits();
	auto fermionTransformation =
			xacc::optionExists("fermion-transformation") ?
					xacc::getOption("fermion-transformation") : "";

	Eigen::VectorXd eigenvalues;
	if (xacc::optionExists("diag-number-symmetry") && 
			xacc::optionExists("n-electrons")) {
		int nElectrons = std::stoi(xacc::getOption("n-electrons"));

		// Generate all n-qubit bitstrings with n-electron
		// bits set
		std::string initBitString = "";
		for (int i = 0; i < nQubits - nElectrons; i++)
			initBitString += "0";
		for (int i = 0; i < nElectrons; i++)
			initBitString += "1";

		// Create occupation basis bit strings
		std::vector<std::string> bitStrings;
		do {
			std::cout << initBitString << "\n";
			bitStrings.push_back(initBitString);
		} while (std::next_permutation(initBitString.begin(),
				initBitString.end()));

		// Transform bit strings from occupation basis
		// to bravyi kitaev basis if needed
		if (fermionTransformation == "bk") {

			// Build up BK transformation matrix
			Eigen::MatrixXi B(1,1); B(0,0) = 1;
			while (true) {

				Eigen::MatrixXi oldB = B;
				Eigen::MatrixXi ones = Eigen::MatrixXi::Ones(1,oldB.cols());
				B.resize(oldB.rows()*2, oldB.cols()*2);
				B.setZero();

				// Set 2nd quadrant
				B.block(0,0,oldB.rows(), oldB.cols()) = oldB;

				// Set 4th quadrant
				B.block(oldB.rows(), oldB.cols(), oldB.rows(), oldB.cols()) = oldB;

				// Set the ones in the first quadrant
				B.block(0, B.cols()-oldB.cols(), ones.rows(), ones.cols()) = ones;

				if (B.rows() == nQubits) {
					break;
				} else if (B.rows() > nQubits) {
					std::cout << "Getting Subblock:\n" << B << "\n\n";
					Eigen::MatrixXi subB = B.block(B.rows()-nQubits, B.cols()-nQubits, nQubits, nQubits);
//					subB.block(nQubits-1, 0, 1, nQubits) = Eigen::MatrixXi::Ones(1,nQubits);
					B = subB;
					break;
				}
			}

			std::cout << "SeeleyBK=\n" << B << "\n\n";

			int end = nQubits-1;
			for (int start = 0; ; start++) {
				B.row(start).swap(B.row(end));
				end--;
				if (start == end) break;
			}

			end = nQubits -1;
			for (int start = 0; ; start++) {
				B.col(start).swap(B.col(end));
				end--;
				if (start == end) break;
			}

			std::cout << "TranterBK:\n" << B << "\n";

			std::vector<std::string> newBitStrings;
			for (auto& bs : bitStrings) {
				Eigen::VectorXi x(nQubits);
				for (int i = 0; i < bs.length(); i++) {
					x(i) = bs[i] == '1' ? 1 : 0;
				}

				Eigen::VectorXi y = B * x;
				for (int i = 0; i < y.size(); i++) y(i) %= 2;

				std::string newbitstring = "";
				for (int i = 0; i < y.size(); i++) {
					newbitstring += y(i) == 0 ? '0' : '1';
				}
				newBitStrings.push_back(newbitstring);
			}

			bitStrings.clear();
			bitStrings = newBitStrings;
		}

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
		std::size_t dim = 1;
		std::size_t two = 2;
		for (int i = 0; i < nQubits; i++)
			dim *= two;
		
		auto getBitStrForIdx = [&](std::uint64_t i) {
			std::stringstream s;
			for (int k = nQubits - 1; k >= 0; k--) s << ((i >> k) & 1);
			return s.str();
		};
		Eigen::MatrixXcd A(dim, dim);
		A.setZero();
		for (std::uint64_t myRow = 0; myRow < dim; myRow++) {
			auto rowBitStr = getBitStrForIdx(myRow);
			auto results = hamiltonian.computeActionOnBra(rowBitStr);
			for (auto& result : results) {
				std::uint64_t k = std::stol(result.first, nullptr, 2);
				A(myRow, k) += result.second;
			}
		}

		es.compute(A);
		eigenvalues = es.eigenvalues();
	}

	gsReal = eigenvalues(0);
	std::stringstream ss;
	ss << std::setprecision(12) << gsReal;
	xacc::info("Ground State Energy of Hamiltonian = " + ss.str());
	return std::real(gsReal);
}

std::pair<double, Eigen::VectorXcd> EigenDiagonalizeBackend::diagonalizeWithGroundState(std::shared_ptr<VQEProgram> prog) {
    auto ground = diagonalize(prog);
    return {ground, es.eigenvectors().col(0)};
}
}
}
