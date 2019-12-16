#pragma once 
#include "eos.hpp"
#include "poly34.h"
#include <algorithm>    // std::sort

constexpr double TwoToOneAndHalf = 2.828427125;
constexpr double TwoToOneHalf = 1.414213562;


double fugacity(double p, double z, double A, double B) {
	double k = log((z + (TwoToOneHalf + 1) * B) / (z - (TwoToOneHalf - 1) * B));
	return p * exp(z - 1 - log(z - B) - (A / (TwoToOneAndHalf * B)) * k);
};

Equilibrium PengRonbinson(Gas& gas, double T, double p) {
	double Pc = gas.getPc();
	double Tc = gas.getTc();
	double omega = gas.getOmega();
	double omegaB = 0.07780;
	double omegaA = 0.45724;
	double R = 10.732;

	double Tr = T / Tc;
	double b = omegaB * (R * Tc) / Pc;
	double ac = omegaA * (R * R * Tc * Tc) / Pc;
	double alfa = pow(1 + (0.37464 + 1.5422 * omega - 0.2699 * omega * omega) * (1 - sqrt(Tr)), 2);
	double at = ac * alfa;

	double A = (at * p) / (R * R * T * T);
	double B = (b * p) / (R * T);

	// Build Cubic Eos z^3 + c0*z^2 +c1*z + c2* = 0 

	double c0 = -(1 - B);
	double c1 = (A - 2 * B - 3 * B * B);
	double c2 = -(A * B - B * B - B * B * B);

	double roots[3];

	auto nroots = SolveP3(roots, c0, c1, c2);
	std::vector<double> solutions(roots, roots + nroots);

	std::sort(solutions.begin(), solutions.end());

	Equilibrium eq;
	

	if (solutions.size() == 3) {
		auto fL = fugacity(p, solutions[0], A, B);
		auto fG = fugacity(p, solutions[2], A, B);
		eq.fugacity.push_back(fL);
		eq.fugacity.push_back(fG);
		eq.zFactors.push_back(solutions[0]);
		eq.zFactors.push_back(solutions[2]);
	}
	else
	{
		auto fL = fugacity(p, solutions[0], A, B);
		eq.fugacity.push_back(fL);
		eq.zFactors.push_back(solutions[0]);
	}
	
	return eq;
}
