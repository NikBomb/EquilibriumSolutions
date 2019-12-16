#include <vector>
#include <functional>
#include "gases.hpp"
#include "eos.hpp" 


/*Example 15.1  -> Non ideal solutions using K factor*/

int main() {

	IsoC4 isoC4{};
	auto EosSol = PengRonbinson(isoC4, 190 + 459.67, 228.79);
};