#pragma once

#include "gases.hpp"


struct Equilibrium {
	std::vector<double> zFactors;
	std::vector<double> fugacity;
};

struct MixtureEquilibrium {
	std::vector<double> liqComp;
	std::vector<double> gasComp;
	double zFactorLiquid;
	double zFactorGas;
	double ngBar;
	double nLBar;
};

MixtureEquilibrium PengRonbinson(GasMixtures& mixture, double T, double p);
Equilibrium PengRonbinson(Gas& gas, double T, double p);