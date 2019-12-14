#include "bisection.hpp"

double bisection(double a, double b, std::function<double(double)> fct)
{
	double e = 1e-6;
	double c;
	if (fct(a) * fct(b) >= 0)
	{
		printf("Incorrect a and b");
		return 0;
	}

	c = a;

	while ((b - a) >= e)
	{
		c = (a + b) / 2;
		if (fct(c) == 0.0) {
			break;
		}
		else if (fct(c) * fct(a) < 0) {
			b = c;
		}
		else {
			a = c;
		}
	}
	return a;
}