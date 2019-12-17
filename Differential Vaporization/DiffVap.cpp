#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

#include "bisection.hpp"
#include "nonIdealSolutions.hpp"
#include "gases.hpp"

/*Example 12.8 -> Differential Vaporization --> result from book n = 0.97 moles*/

int main() {

	// Hypothesis one mole of gas

	C1 c1 = C1();
	C3 c3 = C3();
	nC5 nc5 = nC5();

	double pIn = 600;
	double pFin = 400;
	double pIncr = 0.2;
	double T = 80 + 459.67;

	std::vector<double> z{ 0.160, 0.380, 0.460 };
	std::vector<Gas*> gases{ &c1, &c3, &nc5 };
	std::vector<double> k(3);
	std::vector<double> nLiquidMoles(3);

	double n = 1.0; /*Initial number of moles*/
	double initialbpp = bubblePointPressure(z, gases, T);
	double p = bubblePointPressure(z, gases, T) - pIncr;
	double dN = 0;
	std::vector<double> gasMoles(z.size());
	std::vector<double> newComposition(z.size());


	for (auto i = 0; i < z.size(); i++) {
		nLiquidMoles[i] = z[i] * n;
	}

	while (p > pFin) {
		for (auto i = 0; i < z.size(); i++) {
			k[i] = gases[i]->kEq(p, T);
		}

		std::function<double(double)> bindGasComp = [&z, &k, &p](double nl) { auto gComp = gasComposition(z, k, p, nl);
		return std::accumulate(gComp.begin(), gComp.end(), 0.0) - 1.0;; };



		auto nLBar = (bisection(0.01, 1, bindGasComp));
		auto nL = nLBar * n;
		auto nG = (1 - nLBar) * n;
		n -= nG;
		auto gComp = gasComposition(z, k, p, nLBar);
		newComposition = liquidComposition(z, k, p, 1 - nLBar);

		z = newComposition;
		auto p2 = bubblePointPressure(z, gases, T);
		p -= pIncr;
	}
}