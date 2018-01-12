#include "EfficientJW.hpp"
#include "XACC.hpp"


namespace xacc {
namespace vqe {

std::shared_ptr<IR> EfficientJW::transform(
		std::shared_ptr<IR> ir) {

	boost::mpi::communicator world;

	auto fermiKernel = ir->getKernels()[0];

	result.clear();

	int myStart = 0;
	int myEnd = fermiKernel->nInstructions();

	auto instructions = fermiKernel->getInstructions();
	auto instVec = std::vector<InstPtr>(instructions.begin(), instructions.end());

	auto start = std::clock();
	// Loop over all Fermionic terms...
	for (int z = myStart; z < myEnd; ++z) {

		auto f = instVec[z];

		auto coeff = boost::get<std::complex<double>>(
				f->getParameter(f->nParameters() - 1));

		// Get the creation or annihilation sites
		auto termSites = f->bits();

		if (termSites.size() == 2) {

			int i = termSites[0];
			int j = termSites[1];

			int pmin = std::min(i,j);
			int pmax = std::max(i,j);

			PauliOperator xi({{i,"X"}}, std::complex<double>(0.5, 0)),
					yi({{i,"Y"}}, std::complex<double>(0,.5));
			PauliOperator xj({{j,"X"}}, std::complex<double>(0.5, 0)),
					yj({{j,"Y"}}, std::complex<double>(0,-.5));

			// Construct 1/2 ( xi + iyi) = s_i^+
			xi += yi;

			// Construct 1/2 ( xj - iyj) = s_j^-
			xj += yj;

			std::map<int, std::string> zpm;
			for (int p = pmin; p < pmax; ++p) zpm[p] = "Z";

			result += coeff * xi * PauliOperator(zpm) * xj;

		} else if (termSites.size() == 4) {
			int i = termSites[0];
			int j = termSites[1];
			int k = termSites[2];
			int l = termSites[3];

			int pmin = std::min(j, k);
			int pmax = std::max(j, k);

			PauliOperator xi( { { i, "X" } }, std::complex<double>(0.5, 0)),
					yi( { { i, "Y" } }, std::complex<double>(0, .5));
			PauliOperator xj( { { j, "X" } }, std::complex<double>(0.5, 0)),
					yj( { { j, "Y" } }, std::complex<double>(0, .5));
			PauliOperator xk( { { k, "X" } }, std::complex<double>(0.5, 0)),
					yk( { { k, "Y" } }, std::complex<double>(0, -.5));
			PauliOperator xl( { { l, "X" } }, std::complex<double>(0.5, 0)),
					yl( { { l, "Y" } }, std::complex<double>(0, -.5));

			// Construct 1/2 ( xi + iyi) = s_i^+
			xi += yi;

			// Construct 1/2 ( xj + iyj) = s_j^+
			xj += yj;

			// Construct 1/2 ( xk - iyk) = s_k^+
			xk += yk;

			// Construct 1/2 ( xl - iyl) = s_l^+
			xl += yl;

			// Construct s_k^- * s_l^-
			xk *= xl;

			// Construct s_i^+ * s_j^+
			xi *= xj;

			std::map<int, std::string> zpm;
			for (int p = pmin; p < pmax; ++p) zpm[p] = "Z";

			result += coeff * xi * PauliOperator(zpm) * xk;
		} else if (termSites.size() == 0) {
			result += PauliOperator(coeff);
		}
	}

	std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << "\n";
	return generateIR();
}

}
}

