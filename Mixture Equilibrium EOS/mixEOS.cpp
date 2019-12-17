#include <vector>
#include <functional>
#include "gases.hpp"
#include "eos.hpp" 


/*Example 15.4  -> Non ideal solutions using K factor*/

int main() {
  
	C1 c1{};
	nC4 c4{};
	nC10 c10{};
	std::vector<Gas*> g{ &c1, &c4, &c10 };
	std::vector<double> c{ 0.5301, 0.1055, 0.3644 };
	GasMixtures gm(g, c);
	MixtureEquilibrium q = PengRonbinson(gm, 160 + 459.67, 1000);

};