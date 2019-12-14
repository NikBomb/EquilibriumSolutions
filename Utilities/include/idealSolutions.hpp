#pragma once
#include <vector>
#include <algorithm>
#include <numeric>

std::vector<double> liquidComposition(std::vector<double> z, std::vector<double> pv, double p, double nl)
{
	std::vector<double> partialComp;

	for (int i = 0; i < pv.size(); ++i) {
		partialComp.push_back(z[i] / (1 + nl * (p / pv[i] - 1)));
	}
	return partialComp;
}

std::vector<double> gasComposition(std::vector<double> z, std::vector<double> pv, double p, double ng)
{
	std::vector<double> partialComp;

	for (int i = 0; i < pv.size(); ++i) {
		partialComp.push_back(z[i] / (1 + ng * (pv[i] / p - 1)));
	}
	return partialComp;
}

double dewPoint(std::vector<double> z, std::vector<double> pv) {
	std::vector<double> zp(z.size());
	std::transform(z.begin(), z.end(), pv.begin(), zp.begin(), [](double zj, double pvj) {return zj / pvj; });
	return 1 / std::accumulate(zp.begin(), zp.end(), 0.0);
};

double bubblePoint(std::vector<double> z, std::vector<double> pv) {
	std::vector<double> zp(z.size());
	std::transform(z.begin(), z.end(), pv.begin(), zp.begin(), [](double zj, double pvj) {return zj * pvj; });
	return std::accumulate(zp.begin(), zp.end(), 0.f);
};
