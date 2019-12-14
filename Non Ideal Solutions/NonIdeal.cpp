#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

#include "bisection.hpp"
#include "nonIdealSolutions.hpp"
#include "gases.hpp"

/*Example 12.4 -> Non ideal solutions using K factor*/

int main() {

	std::vector<double> z{ 0.610, 0.280, 0.110 };

	C3 c3{};
	nC4 nc4{};
	nC5 nc5{};

	//std::vector<double> k{ 1.55, 0.592, 0.236 };
	double p = 200;
	double T = 150 + 459.67;

	//std::vector<double> k{1.55, 0.592, 0.236};
	std::vector<double> k;
	k.push_back(c3.kEq(p, T));
	k.push_back(nc4.kEq(p, T));
	k.push_back(nc5.kEq(p, T));

	std::vector<Gas*> gases{ &c3, &nc4, &nc5 };

	std::function<double(double)> bindGasComp = [&z, &k, &p](double nl) { auto lComp = gasComposition(z, k, p, nl);
	return std::accumulate(lComp.begin(), lComp.end(), 0.0) - 1.0; };
	auto nL = bisection(0.01, 1, bindGasComp);
	auto liquidCompositions = liquidComposition(z, k, p, nL);
	auto gComp = gasComposition(z, k, p, 1 - nL);

	auto bbp = bubblePointPressure(z, gases, T);
	auto dpp = dewPointPressure(z, gases, T);

}