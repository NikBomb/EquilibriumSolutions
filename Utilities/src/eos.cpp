#pragma once 

#include <algorithm>    // std::sort
#include <vector>
#include <math.h>
#include "eos.hpp"
#include "poly34.h"
#include "nonIdealSolutions.hpp"

#include <typeindex>
#include <typeinfo>


constexpr double TwoToOneAndHalf = 2.828427125;
constexpr double TwoToOneHalf = 1.414213562;

struct PengRobinsonData {
	double omega;
	double omegaB;
	double omegaA;
	double R;

	double b;
	double ac;
	double alfa;
	double at;
};

PengRobinsonData InitializePr(Gas* gas, double T) {
	PengRobinsonData Pr;

	double Pc = gas->getPc();
	double Tc = gas->getTc();
	Pr.omega = gas->getOmega();
	Pr.omegaB = 0.07780;
	Pr.omegaA = 0.45724;
	Pr.R = 10.732;

	double Tr = T / Tc;
	Pr.b = Pr.omegaB * (Pr.R * Tc) / Pc;
	Pr.ac = Pr.omegaA * (Pr.R * Pr.R * Tc * Tc) / Pc;
	Pr.alfa = pow(1 + (0.37464 + 1.5422 * Pr.omega - 0.2699 * Pr.omega * Pr.omega) * (1 - sqrt(Tr)), 2);
	Pr.at = Pr.ac * Pr.alfa;

	return Pr;
}

double fugacity(double p, double z, double A, double B) {
	double k = log((z + (TwoToOneHalf + 1) * B) / (z - (TwoToOneHalf - 1) * B));
	return p * exp(z - 1 - log(z - B) - (A / (TwoToOneAndHalf * B)) * k);
};


double computeZfactor(std::vector<double>& comp, GasMixtures& mixture, std::vector<PengRobinsonData>& Pr, std::vector<double> k, double p, double T) {
	double b = 0;
	double atG = 0;

	for (int i = 0; i < comp.size(); ++i) {
		b += comp[i] * Pr[i].b;
		for (int j = 0; j < k.size(); ++j) {
			double interaction = mixture.getBinaryInteraction(mixture.gases[i], mixture.gases[j]);
			atG += comp[i] * comp[j] * sqrt(Pr[i].at * Pr[j].at) * (1 - interaction);
		}
	}

	double A = atG * p / (Pr[0].R * Pr[0].R * T * T);
	double B = b * p / (Pr[0].R * T);

	double c0 = -(1 - B);
	double c1 = (A - 2 * B - 3 * B * B);
	double c2 = -(A * B - B * B - B * B * B);

	double roots[3];

	auto nroots = SolveP3(roots, c0, c1, c2);
	
	return roots[0];

}
Equilibrium PengRonbinson(GasMixtures& mixture, double T, double p) {

	std::vector<PengRobinsonData> Pr;
	std::vector<double> k;
	std::vector<double> z = mixture.comp;


	std::transform(mixture.gases.begin(), mixture.gases.end(), std::back_inserter(Pr), [&T](auto g) {return InitializePr(g, T); });
	std::transform(mixture.gases.begin(), mixture.gases.end(), std::back_inserter(k), [&T, &p](auto g) {return g->kEq(p, T); });


	std::function<double(double)> bindLiquidComp = [&z, &k, &p, &T](double ng) { auto lComp = liquidComposition( z, k, p, ng);
	return std::accumulate(lComp.begin(), lComp.end(), 0.0) - 1.0; };

	auto ng = bisection(0.01, 1, bindLiquidComp);
	std::vector<double> gComp = gasComposition(z, k, p, 1 - ng);
	std::vector<double> lComp = liquidComposition(z, k, p, ng);

	double zG = computeZfactor(gComp, mixture, Pr, k, p, T);
	double zL = computeZfactor(lComp, mixture, Pr, k, p, T);
	
	return Equilibrium{};
}

Equilibrium PengRonbinson(Gas& gas, double T, double p) {

	auto Pr = InitializePr(&gas, T);

	double A = (Pr.at * p) / (Pr.R * Pr.R * T * T);
	double B = (Pr.b * p) / (Pr.R * T);

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

