#include "EfficientJW.hpp"
#include "XACC.hpp"


namespace xacc {
namespace vqe {

std::shared_ptr<IR> EfficientJW::transform(
		std::shared_ptr<IR> ir) {

	std::complex<double> imag(0,1);
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

		auto coeff =
				f->getParameter(f->nParameters() - 2).as<std::complex<double>>();

		auto fermionVar = f->getParameter(f->nParameters() - 1).as<std::string>();

		// Get the creation or annihilation sites
		auto termSites = f->bits();

		if (termSites.size() == 2) {

			int i = termSites[0];
			int j = termSites[1];

			int pmin = std::min(i,j);
			int pmax = std::max(i,j);

			PauliOperator Sxi({{i,"X"}}, 0.5), Syi({{i,"Y"}}, 0.5);
			PauliOperator Sxj({{j,"X"}}, 0.5), Syj({{j,"Y"}}, 0.5);

			PauliOperator sPlusI = Sxi - imag * Syi;
			PauliOperator sMinusJ = Sxj + imag * Syj;

			std::map<int, std::string> zpm;
			int parity = 1;
			for (int p = pmin; p < pmax-1; ++p) {
				zpm[p] = "Z";
				parity *= -1;
			}

			result += coeff * sPlusI * PauliOperator(zpm, parity) * sMinusJ;

		} else if (termSites.size() == 4) {
			int i = termSites[0];
			int j = termSites[1];
			int k = termSites[2];
			int l = termSites[3];

			int pmin = std::min(j, k);
			int pmax = std::max(j, k);

			PauliOperator Sxi({{i,"X"}}, 0.5), Syi({{i,"Y"}}, 0.5);
			PauliOperator Sxj({{j,"X"}}, 0.5), Syj({{j,"Y"}}, 0.5);
			PauliOperator Sxk({{k,"X"}}, 0.5), Syk({{k,"Y"}}, 0.5);
			PauliOperator Sxl({{l,"X"}}, 0.5), Syl({{l,"Y"}}, 0.5);

			PauliOperator sPlusI = Sxi - imag * Syi;
			PauliOperator sPlusJ = Sxj - imag * Syj;
			PauliOperator sMinusK = Sxk + imag * Syk;
			PauliOperator sMinusL = Sxl + imag * Syl;

			std::map<int, std::string> zpm;
			int parity = 1;
			for (int p = pmin; p < pmax-1; ++p) {
				zpm[p] = "Z";
				parity *= -1;
			}

			result += coeff * sPlusI * sPlusJ * PauliOperator(zpm, parity) * sMinusK * sMinusL;
		} else if (termSites.size() == 0) {
			result += PauliOperator(coeff);
		}
	}

	std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << "\n";
	return result.toXACCIR();
}

}
}

