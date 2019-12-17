#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

#include "gases.hpp"
#include "eos.hpp"

/*Example 12.8 -> Differential Vaporization --> result from book n = 0.97 moles*/

int main() {

	// Hypothesis one mole of gas

	C1 c1 = C1();
	C3 c3 = C3();
	nC5 nc5 = nC5();

	double pIn = 400;
	double pFin = 400;
	double pIncr = 5;
	double T = 80 + 459.67;

	std::vector<double> z{ 0.160, 0.380, 0.460 };
	std::vector<Gas*> gases{ &c1, &c3, &nc5 };
	std::vector<double> k(3);	
	GasMixtures gm(gases, z);
	

	double n = 1.0; /*Initial number of moles*/
	double initialbpp = 585.45;
	double p = initialbpp - pIncr;
	double dN = 0;
	std::vector<double> gasMoles(z.size());
	std::vector<double> newComposition(z.size());




	while (p > pFin) {
		
		MixtureEquilibrium q = PengRonbinson(gm, T , p);

		auto nLBar = q.nLBar;
		auto nL = nLBar * n;
		auto nG = (1 - nLBar) * n;
		n -= nG;
		auto gComp =  q.gasComp;
		newComposition = q.liqComp;

		gm.comp = newComposition;
		
		p -= pIncr;
	}
}