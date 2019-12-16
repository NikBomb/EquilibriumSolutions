#pragma once

#include "gases.hpp"


struct Equilibrium {
	std::vector<double> zFactors;
	std::vector<double> fugacity;
};

Equilibrium PengRonbinson(GasMixtures& mixture, double T, double p);
Equilibrium PengRonbinson(Gas& gas, double T, double p);