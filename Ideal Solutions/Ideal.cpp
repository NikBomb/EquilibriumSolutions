#include<stdio.h>
#include <functional>
#include <vector>
#include <algorithm>
#include <numeric>
#include <assert.h>
#include "bisection.hpp"
#include "idealSolutions.hpp"

//Example 12.1 - 12.3: Flash Calculation
/*
Calculate composition and quantities of the gas and liquid when 1 mole of the given mixture is brought to equilibrium.
*/
int main()
{
	std::vector<double> z{ 0.610, 0.280, 0.110 };
	std::vector<double> pv{ 350, 105, 37 };
	double p = 200;

	std::function<double(double)> bindliqComp = [&z, &pv, &p](double nl) { auto lComp = liquidComposition(z, pv, p, nl);
	return std::accumulate(lComp.begin(), lComp.end(), 0.0) - 1.0;; };
	auto nL = bisection(0.01, 1, bindliqComp);
	auto liquidCompositions = liquidComposition(z, pv, p, nL);
	auto gComp = gasComposition(z, pv, p, 1 - nL);
	auto pb = bubblePoint(z, pv);
	auto sp = dewPoint(z, pv);
}