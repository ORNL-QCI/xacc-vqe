#include "SlepcGroundStateEnergyCalculator.hpp"
#include <petscmat.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <slepceps.h>
#include <unsupported/Eigen/KroneckerProduct>
#include "SpinInstruction.hpp"

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <memory>

using IndexPair = std::pair<std::uint64_t, std::uint64_t>;

namespace xacc {
namespace vqe {
double SlepcGroundStateEnergyCalculator::computeGroundStateEnergy(
		CompositeSpinInstruction& inst, const int nQubits) {
	boost::mpi::communicator world;
	int rank = world.rank();
	int nRanks = world.size();
	double gsReal;
	static char help[] = "";
	std::vector<std::string> argvVec;
	std::vector<char*> cstrs;
	argvVec.insert(argvVec.begin(), "appExec");
	for (auto& s : argvVec) {
		cstrs.push_back(&s.front());
	}

	int argc = argvVec.size();
	auto argv = cstrs.data();

	std::size_t dim = 1;
	std::size_t two = 2;
	for (int i = 0; i < nQubits; i++)
		dim *= two;

	std::unordered_map<IndexPair, std::complex<double>, boost::hash<IndexPair> > nonZeros;

	// We know there will be diagonal enties
	for (int i = 0; i < dim; i++)
		nonZeros.insert(
				std::make_pair(IndexPair { i, i }, std::complex<double>(0, 0)));

	// Generate all bit strings
	std::vector<std::string> bitStrings;
	for (std::uint64_t j = 0; j < dim; j++) {

		std::stringstream s;
		for (int k = nQubits - 1; k >= 0; k--) {
			s << ((j >> k) & 1);
		}
		bitStrings.push_back(s.str());
	}

	if (rank == 0) XACCInfo("Building up list of nonZeros.");
	for (int i = 0; i < inst.nInstructions(); i++) {
		std::shared_ptr<SpinInstruction> spinInst = std::dynamic_pointer_cast<
				SpinInstruction>(inst.getInstruction(i));

		IndexPair p;
		std::pair<std::string, std::complex<double>> newBitStrCoeff;

		if (spinInst->isIdentity()) {
			for (auto& kv : nonZeros) {
				if (kv.first.first == kv.first.second) {
					kv.second += spinInst->coefficient;
				}
			}
		} else if (spinInst->isDiagonal()) {
			for (std::uint64_t j = 0; j < dim; j++) {
				p = std::make_pair(j, j);
				newBitStrCoeff = spinInst->computeActionOnBits(bitStrings[j]);
				nonZeros[p] += newBitStrCoeff.second;
			}
		} else {
			for (std::uint64_t j = 1; j < dim; j++) {
				newBitStrCoeff = spinInst->computeActionOnBits(bitStrings[j]);
				std::uint64_t i = std::stol(newBitStrCoeff.first, nullptr, 2);
				p = std::make_pair(i, j);
				if (nonZeros.find(p) != nonZeros.end()) {
					nonZeros[p] += newBitStrCoeff.second;
				} else {
					nonZeros.insert(std::make_pair(p, newBitStrCoeff.second));
				}
			}
		}
	}
	if (rank == 0) XACCInfo("Done building up list of nonZeros.");

	Mat A;
	EPS eps;
	EPSType type;
	PetscInt n = dim, Istart, Iend;

	SlepcInitialize(&argc, &argv, (char*) 0, help);

	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
	MatSetFromOptions(A);
	MatSetUp(A);

	MatGetOwnershipRange(A, &Istart, &Iend);
	for (auto& kv : nonZeros) {
		int r = kv.first.first;
		int col = kv.first.second;
		double val = std::real(kv.second);
		if (std::fabs(val) > 1e-12) {
			// Only add our rank's values
			if (r >= Istart && r < Iend) {
				MatSetValue(A, r, col, val, INSERT_VALUES);
			}
		}
	}

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	MatView(A, PETSC_VIEWER_STDOUT_WORLD);
	EPSCreate(PETSC_COMM_WORLD, &eps);
	EPSSetOperators(eps, A, NULL);
	EPSSetProblemType(eps, EPS_HEP);
	EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
	EPSSetType(eps, EPSLANCZOS);
	if (rank == 0)
		XACCInfo("Starting SLEPc EPS Solve");
	EPSSolve(eps);
	if (rank == 0)
		XACCInfo("Done with SLEPc EPS Solve");
	EPSGetEigenpair(eps, 0, &gsReal, NULL, NULL, NULL);
	EPSDestroy(&eps);
	MatDestroy(&A);
	SlepcFinalize();

	return std::real(gsReal);
}

}
}

//	using SparseComplexMatrix = Eigen::SparseMatrix<std::complex<double>>;
//	SparseComplexMatrix spMat = inst.toSparseMatrix(nQubits).pruned();
//	spMat.makeCompressed();
//
//	using SparseMat = Eigen::SparseMatrix<std::complex<double>>;
//	using Triplet = Eigen::Triplet<std::complex<double>>;
//
//	std::complex<double> i(0, 1);
//	SparseMat z(2, 2), x(2, 2), y(2, 2), I(2, 2);
//	std::vector<Triplet> zCoeffs{Triplet(0,0,1), Triplet(1,1,-1)};
//	std::vector<Triplet> xCoeffs{Triplet(0,1,1), Triplet(1,0,1)};
//	std::vector<Triplet> yCoeffs{Triplet(0,1,-i), Triplet(1, 0, i)};
//
//	z.setFromTriplets(zCoeffs.begin(), zCoeffs.end());
//	x.setFromTriplets(xCoeffs.begin(), xCoeffs.end());
//	y.setFromTriplets(yCoeffs.begin(), yCoeffs.end());
//
//	std::size_t dim = 1;
//	std::size_t two = 2;
//	for (int i = 0; i < nQubits; i++)
//		dim *= two;
//
//	std::vector<SparseMat> Zs, Xs, Ys;
//
//	if (rank == 0) XACCInfo("Building Zs, Xs, Ys.");
//
//	for (int i = 0; i < nQubits; i++) {
//
//		std::size_t d1 = 1, d2 = 1;
//		for (int j = 0; j < i; j++) d1 *= two;
//		for (int j = 0; j < nQubits-i-1; j++) d2 *= two;
//
//		if (rank == 0) XACCInfo("Creating Matrices for qbit " + std::to_string(i));
//
//		SparseMat i1(d1,d1), iNi(d2,d2);
//		i1.setIdentity();
//		iNi.setIdentity();
//
//		SparseMat zi =
//				Eigen::kroneckerProduct(i1,
//						Eigen::kroneckerProduct(z, iNi).eval().pruned()).eval().pruned();
//		zi.makeCompressed();
//		Zs.push_back(zi);
//
//		SparseMat xi =
//				Eigen::kroneckerProduct(i1,
//						Eigen::kroneckerProduct(x, iNi).eval().pruned()).eval().pruned();
//		xi.makeCompressed();
//		Xs.push_back(xi);
//
//		SparseMat yi =
//				Eigen::kroneckerProduct(i1,
//						Eigen::kroneckerProduct(y, iNi).eval().pruned()).eval().pruned();
//		yi.makeCompressed();
//		Ys.push_back(yi);
//	}
//	if (rank == 0) XACCInfo("Done building Zs, Xs, Ys.");
//
//	if (rank == 0) XACCInfo("Constructing Hamiltonian matrix with " + std::to_string(inst.nInstructions()) + " terms.");
//	Eigen::SparseMatrix<std::complex<double>> ham (dim, dim);
//	for (int i = 0; i < inst.nInstructions(); i++) {
//		auto j = std::dynamic_pointer_cast<SpinInstruction>(inst.getInstruction(i));
//		std::vector<std::pair<int, std::string>> terms = j->getTerms();
//
//		if (terms.size() == 1 && terms[0] == std::pair<int, std::string> { 0,
//				"I" }) {
//			SparseMat id(dim,dim);
//			id.setIdentity();
//			ham += j->coefficient * id;
//			continue;
//		}
//
//		SparseMat local(dim, dim);
//		local.setIdentity();
//		for (int k = 0; k < terms.size(); k++) {
//			auto t = terms[k];
//			auto qbit = t.first;
//			auto gate = t.second;
//
//			SparseMat tmp;
//			if (gate == "Z") {
//				tmp = Zs[qbit];
//			} else if (gate == "Y") {
//				tmp = Ys[qbit];
//			} else if (gate == "X") {
//				tmp = Xs[qbit];
//			} else {
//				XACCError("Invalid gate name - " + gate);
//			}
//
//			local = local*tmp;
//			local.makeCompressed();
//		}
//
//		ham += j->coefficient * local;
//	}
//	if (rank == 0) XACCInfo("Done constructing Hamiltonian matrix with " + std::to_string(inst.nInstructions()) + " terms.");
//	if (rank == 0) XACCInfo("There are " + std::to_string(ham.nonZeros()) + " non zero elements.");
//
//	if (rank == 0) std::cout << ham << "\n";
