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

struct MixtureConstants {
	double zFactor;
	double at;
	double b;
	double A;
	double B;
	std::vector<double> composition;
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

double fugacityCoefficientMixtures(double z, double A, double B, double Aj, double Bj) {
	double k = log((z + (TwoToOneHalf + 1) * B) / (z - (TwoToOneHalf - 1) * B));
	return exp(-log(z - B) + (z - 1) * Bj - (A / (TwoToOneAndHalf * B)) * (Aj - Bj) * k);
};

MixtureConstants computeMixtureConstants(std::vector<double>& comp, GasMixtures& mixture, std::vector<PengRobinsonData>& Pr, std::vector<double> k, double p, double T) {
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

	MixtureConstants mc;
	mc.zFactor = roots[0];
	mc.at = atG;
	mc.b = b;
	mc.A = A;
	mc.B = B;

	return mc;

}
MixtureEquilibrium PengRonbinson(GasMixtures& mixture, double T, double p) {

	std::vector<PengRobinsonData> Pr;
	std::vector<double> kt;
	std::vector<double> kc(mixture.gases.size());
	double err = 100;
	double tol = 1e-10;
	double ng;


	std::vector<double> z = mixture.comp;
	MixtureConstants mc[2];

	std::transform(mixture.gases.begin(), mixture.gases.end(), std::back_inserter(Pr), [&T](auto g) {return InitializePr(g, T); });
	std::transform(mixture.gases.begin(), mixture.gases.end(), std::back_inserter(kt), [&T, &p](auto g) {return g->kEq(p, T); });

	while (err > tol) {
		std::function<double(double)> bindgasComp = [&z, &kt, &p, &T](double nL) { auto lComp = gasComposition(z, kt, p, nL);
		return std::accumulate(lComp.begin(), lComp.end(), 0.0) - 1.0; };

		auto nL = bisection(0.01, 1, bindgasComp);
		ng = 1 - nL;
		std::vector<double> gComp = gasComposition(z, kt, p, 1 - ng);
		std::vector<double> lComp = liquidComposition(z, kt, p, ng);

		mc[0] = computeMixtureConstants(gComp, mixture, Pr, kt, p, T);
		mc[1] = computeMixtureConstants(lComp, mixture, Pr, kt, p, T);

		mc[0].composition = gComp;
		mc[1].composition = lComp;

		std::vector<std::vector<double>> A;
		std::vector<std::vector<double>> B;
		std::vector<std::vector<double>> fug;


		for (int f = 0; f < 2; ++f) {
			double at = mc[f].at;
			double b = mc[f].b;
			std::vector<double> c = mc[f].composition;
			std::vector<double> Aij;
			std::vector<double> Bij;

			for (int j = 0; j < mixture.gases.size(); ++j) {
				double Aj = (2 * sqrt(Pr[j].at)) / at;
				double Bj = Pr[j].b / mc[f].b;
				double AjSum = 0;
				for (int i = 0; i < mixture.gases.size(); ++i) {
					AjSum += c[i] * sqrt(Pr[i].at) * (1 - mixture.getBinaryInteraction(mixture.gases[i], mixture.gases[j]));
				}
				Aij.push_back(Aj * AjSum);
				Bij.push_back(Bj);
			}
			A.push_back(Aij);
			B.push_back(Bij);
		}

		for (int f = 0; f < 2; ++f) {
			std::vector<double> ff;
			for (int j = 0; j < mixture.gases.size(); ++j) {
				ff.push_back(fugacityCoefficientMixtures(mc[f].zFactor, mc[f].A, mc[f].B, A[f][j], B[f][j]));
			}
			fug.push_back(ff);
		}


		for (int j = 0; j < mixture.gases.size(); ++j) {
			kc[j] = (fug[1][j] / (fug[0][j]));
			err = (kc[j] - kt[j]) * (kc[j] - kt[j]) / (kc[j] * kt[j]);
		}
		kt = kc;
	}

	MixtureEquilibrium mE;

	mE.gasComp = mc[0].composition;
	mE.liqComp = mc[1].composition;
	mE.zFactorGas = mc[0].zFactor;
	mE.zFactorLiquid = mc[1].zFactor;
	mE.ngBar = ng;
	mE.nLBar = 1 - ng;


	return mE;
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

