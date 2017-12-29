#include "SpinInstruction.hpp"


using Triplet = Eigen::Triplet<std::complex<double>>;
using SparseMat = Eigen::SparseMatrix<std::complex<double>>;

xacc::vqe::SpinInstruction operator*(double const& scalar,
		xacc::vqe::SpinInstruction rhs) {
	return rhs *= scalar;
}

xacc::vqe::SpinInstruction operator*(std::complex<double> const& scalar,
		xacc::vqe::SpinInstruction rhs) {
	return rhs *= scalar;
}

#pragma omp declare reduction( * : Eigen::SparseMatrix<std::complex<double>> : omp_out = omp_out * omp_in ) \
  initializer( omp_priv = omp_orig )

namespace xacc {
namespace vqe {
const std::pair<std::string, std::complex<double>> SpinInstruction::computeActionOnKet(const std::string& bitString) {
	std::complex<double> i(0,1);
	std::complex<double> newCoeff = coefficient;
	std::string newBits = bitString;
	for (auto& t : terms) {
		auto idx = t.first;
		auto gate = t.second;
		if (gate == "Z") {
			newCoeff *= newBits[idx] == '1' ? -1 : 1;
		} else if (gate == "X") {
			newBits[idx] = (newBits[idx] == '1' ? '0' : '1');
		} else if (gate == "Y") {
			newCoeff *= newBits[idx] == '1' ? -i : i;
			newBits[idx] = (newBits[idx] == '1' ? '0' : '1');
		}
	}
	return std::make_pair(newBits, newCoeff);
}

const std::pair<std::string, std::complex<double>> SpinInstruction::computeActionOnBra(const std::string& bitString) {
	std::complex<double> i(0,1);
	std::complex<double> newCoeff = coefficient;
	std::string newBits = bitString;

	for (auto& t : terms) {
		auto idx = t.first;
		auto gate = t.second;
		if (gate == "Z") {
			newCoeff *= newBits[idx] == '1' ? -1 : 1;
		} else if (gate == "X") {
			newBits[idx] = (newBits[idx] == '1' ? '0' : '1');
		} else if (gate == "Y") {
			newCoeff *= newBits[idx] == '1' ? i : -i;
			newBits[idx] = (newBits[idx] == '1' ? '0' : '1');
		}
	}
	return std::make_pair(newBits, newCoeff);

}
Eigen::SparseMatrix<std::complex<double>> SpinInstruction::toSparseMatrix(
		const int nQubits) {
	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	SparseMat ham(dim, dim);
	ham.setIdentity();

	if (isIdentity()) {
		return coefficient * ham;
	}

	std::complex<double> i(0, 1);
	SparseMat z(2, 2), x(2, 2), y(2, 2), I(2, 2);
	std::vector<Triplet> zCoeffs { Triplet(0, 0, 1), Triplet(1, 1, -1) };
	std::vector<Triplet> xCoeffs { Triplet(0, 1, 1), Triplet(1, 0, 1) };
	std::vector<Triplet> yCoeffs { Triplet(0, 1, -i), Triplet(1, 0, i) };

	z.setFromTriplets(zCoeffs.begin(), zCoeffs.end());
	x.setFromTriplets(xCoeffs.begin(), xCoeffs.end());
	y.setFromTriplets(yCoeffs.begin(), yCoeffs.end());

	for (auto& kv : terms) {

		auto qbit = kv.first;
		auto gate = kv.second;

		SparseMat tmp;
		if (gate == "Z") {
			tmp = z;
		} else if (gate == "Y") {
			tmp = y;
		} else if (gate == "X") {
			tmp = x;
		} else {
			XACCError("Invalid gate name - " + gate);
		}

		std::size_t d1 = 1, d2 = 1;
		for (int i = 0; i < qbit; i++)
			d1 *= two;
		for (int i = 0; i < nQubits - qbit - 1; i++)
			d2 *= two;

		SparseMat i1(d1, d1), iNi(d2, d2);
		i1.setIdentity();
		iNi.setIdentity();

		SparseMat localU =
				Eigen::kroneckerProduct(i1,
						Eigen::kroneckerProduct(tmp, iNi).eval().pruned()).eval().pruned();
		localU.makeCompressed();

		ham = ham * localU;
	}

	return coefficient * ham;
}

std::vector<int> SpinInstruction::toBinaryVector(const int nQubits) const {

	std::vector<int> bv(2 * nQubits);

	for (auto& t : terms) {
		if (t.second == "X") {

			bv[t.first] = 1;

		} else if (t.second == "Z") {

			bv[nQubits + t.first] = 1;

		} else if (t.second == "Y") {

			bv[t.first] = 1;
			bv[nQubits + t.first] = 1;

		} else if (t.second == "I") {

			return bv;
		}
	}

	return bv;
}

void SpinInstruction::fromBinaryVector(std::vector<int> vec,
		std::complex<double> coeff) {

	coefficient = coeff;
	int nQubits = vec.size() / 2;

	terms.clear();

	for (int i = 0; i < nQubits; i++) {
		if (vec[i] && vec[i + nQubits]) {
			terms.insert( { i, "Y" });
		} else if (vec[i]) {
			terms.insert( { i, "X" });
		}
	}

	for (int i = nQubits; i < vec.size(); i++)
		if (vec[i] && !vec[i - nQubits])
			terms.insert( { i - nQubits, "Z" });

	return;

}

Eigen::SparseMatrix<double> SpinInstruction::toSparseRealMatrix(
		const int nQubits) {
	Eigen::SparseMatrix<std::complex<double>> mat = toSparseMatrix(nQubits);
	std::vector<Eigen::Triplet<double>> triplets;
	for (int k = 0; k < mat.outerSize(); ++k) {
		for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(mat,
				k); it; ++it) {
			if (it.value() != std::complex<double>(0.0, 0.0)) {
				triplets.push_back(
						Eigen::Triplet<double>(it.row(), it.col(),
								std::real(it.value())));
			}
		}
	}
	Eigen::SparseMatrix<double> spRealMat(mat.rows(), mat.cols());
	spRealMat.setFromTriplets(triplets.begin(), triplets.end());
	return spRealMat;
}

Eigen::MatrixXcd SpinInstruction::toMatrix(const int nQubits) {

	int dim = std::pow(2, nQubits);

	Eigen::MatrixXcd ham = Eigen::MatrixXcd::Identity(dim, dim);

	if (isIdentity()) {
		return coefficient * Eigen::MatrixXcd::Identity(dim, dim);
	}

	std::complex<double> i(0, 1);
	Eigen::MatrixXcd z(2, 2), x(2, 2), y(2, 2), I(2, 2);
	z << 1, 0, 0, -1;
	x << 0, 1, 1, 0;
	y << 0, -i, i, 0;
	I << 1, 0, 0, 1;

	for (auto t : terms) {

		std::vector<Eigen::MatrixXcd> productList;
		for (int j = 0; j < nQubits; j++) {
			productList.push_back(I);
		}

		auto qbit = t.first;
		auto gate = t.second;

		Eigen::MatrixXcd tmp;
		if (gate == "Z") {
			tmp = z;
		} else if (gate == "Y") {
			tmp = y;
		} else if (gate == "X") {
			tmp = x;
		} else {
			XACCError("Invalid gate name - " + gate);
		}

		productList.at(qbit) = tmp;

		Eigen::MatrixXcd localU = productList.at(0);
		for (int j = 1; j < productList.size(); j++) {
			localU = Eigen::kroneckerProduct(localU, productList.at(j)).eval();
		}

		ham = ham * localU;
	}

	return coefficient * ham;
}

}
}
