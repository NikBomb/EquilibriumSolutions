#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include "gases.hpp"
#include "bisection.hpp"


std::vector<double> gasComposition(std::vector<double> z, std::vector<double> k, double p, double nl)
{
	std::vector<double> partialComp;

	for (int i = 0; i < k.size(); ++i) {
		partialComp.push_back(z[i] / (1 + nl * (1 / k[i] - 1)));
	}
	return partialComp;
}

std::vector<double> liquidComposition(std::vector<double> z, std::vector<double> k, double p, double ng)
{
	std::vector<double> partialComp;

	for (int i = 0; i < k.size(); ++i) {
		partialComp.push_back(z[i] / (1 + ng * (k[i] - 1)));
	}
	return partialComp;
}

double bubblePointPressure(std::vector<double> z, std::vector<Gas*> gases, double T) {
	double pMin = 100;
	double pMax = 1000;

	std::vector<double> k(gases.size());

	auto bubblePointEquation = [&](double p) {
		std::transform(gases.begin(), gases.end(), k.begin(), [&](auto g) {return g->kEq(p, T); });
		std::vector<double> kTimesz(gases.size());
		std::transform(k.begin(), k.end(), z.begin(), kTimesz.begin(), [&](auto kj, auto zj) {return zj * kj; });
		return std::accumulate(kTimesz.begin(), kTimesz.end(), 0.0) - 1;
	};

	return bisection(pMin, pMax, bubblePointEquation);
};

double dewPointPressure(std::vector<double> z, std::vector<Gas*> gases, double T) {
	double pMin = 30;
	double pMax = 1000;

	std::vector<double> k(gases.size());

	auto bubblePointEquation = [&](double p) {
		std::transform(gases.begin(), gases.end(), k.begin(), [&](auto g) {return g->kEq(p, T); });
		std::vector<double> kTimesz(gases.size());
		std::transform(k.begin(), k.end(), z.begin(), kTimesz.begin(), [&](auto kj, auto zj) {return zj / kj; });
		return std::accumulate(kTimesz.begin(), kTimesz.end(), 0.0) - 1;
	};

	return bisection(pMin, pMax, bubblePointEquation);
};
