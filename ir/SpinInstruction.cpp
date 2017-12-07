#include "SpinInstruction.hpp"

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


namespace xacc{
namespace vqe {

Eigen::SparseMatrix<std::complex<double>> SpinInstruction::toSparseMatrix(const int nQubits) {
	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	SparseMat ham(dim,dim);
	ham.setIdentity();

	if (terms.size() == 1 && terms[0] == std::pair<int, std::string> { 0,
			"I" }) {
		return coefficient
				* ham;
	}

	std::complex<double> i(0, 1);
	SparseMat z(2, 2), x(2, 2), y(2, 2), I(2, 2);
	std::vector<Triplet> zCoeffs{Triplet(0,0,1), Triplet(1,1,-1)};
	std::vector<Triplet> xCoeffs{Triplet(0,1,1), Triplet(1,0,1)};
	std::vector<Triplet> yCoeffs{Triplet(0,1,-i), Triplet(1, 0, i)};

	z.setFromTriplets(zCoeffs.begin(), zCoeffs.end());
	x.setFromTriplets(xCoeffs.begin(), xCoeffs.end());
	y.setFromTriplets(yCoeffs.begin(), yCoeffs.end());

#pragma omp parallel for shared(terms) reduction (*:ham)
	for (int i = 0; i < terms.size(); i++) {

		auto t = terms[i];

		auto qbit = t.first;
		auto gate = t.second;

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
		for (int i = 0; i < qbit; i++) d1 *= two;
		for (int i = 0; i < nQubits-qbit-1; i++) d2 *= two;


		SparseMat i1(d1,d1), iNi(d2,d2);
		i1.setIdentity();
		iNi.setIdentity();

		SparseMat localU =
				Eigen::kroneckerProduct(i1,
						Eigen::kroneckerProduct(tmp, iNi).eval().pruned()).eval().pruned();
		localU.makeCompressed();

		ham = ham * localU;
	}

	return coefficient * ham;

//	std::size_t dim = 1;
//	std::size_t two = 2;
//	for (int i = 0; i < nQubits; i++)
//		dim *= two;
//
//	SparseMat ham(dim,dim);
//	ham.setIdentity();
//
//	if (terms.size() == 1 && terms[0] == std::pair<int, std::string> { 0,
//			"I" }) {
//		return coefficient
//				* ham;
//	}
//
//	std::complex<double> i(0, 1);
//	SparseMat z(2, 2), x(2, 2), y(2, 2), I(2, 2);
//	std::vector<Triplet> zCoeffs{Triplet(0,0,1), Triplet(1,1,-1)};
//	std::vector<Triplet> xCoeffs{Triplet(0,1,1), Triplet(1,0,1)};
//	std::vector<Triplet> yCoeffs{Triplet(0,1,-i), Triplet(1, 0, i)};
//	std::vector<Triplet> ICoeffs{Triplet(0, 0, 1), Triplet(1, 1, 1)};
//
//	z.setFromTriplets(zCoeffs.begin(), zCoeffs.end());
//	x.setFromTriplets(xCoeffs.begin(), xCoeffs.end());
//	y.setFromTriplets(yCoeffs.begin(), yCoeffs.end());
//	I.setFromTriplets(ICoeffs.begin(), ICoeffs.end());
//
//	for (auto t : terms) {
//
//		std::vector<SparseMat> productList;
//		for (int j = 0; j < nQubits; j++) {
//			productList.push_back(I);
//		}
//
//		auto qbit = t.first;
//		auto gate = t.second;
//
//		SparseMat tmp;
//		if (gate == "Z") {
//			tmp = z;
//		} else if (gate == "Y") {
//			tmp = y;
//		} else if (gate == "X") {
//			tmp = x;
//		} else {
//			XACCError("Invalid gate name - " + gate);
//		}
//
//		productList.at(qbit) = tmp;
//
//		SparseMat localU = productList.at(0);
//		for (int j = 1; j < productList.size(); j++) {
//			localU =
//					Eigen::kroneckerProduct(localU, productList.at(j)).eval();
//		}
//
//		ham = ham * localU;
//	}
//
//	return coefficient * ham;

}

}
}
