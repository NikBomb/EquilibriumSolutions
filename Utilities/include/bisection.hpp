#include <functional>

double bisection(double a, double b, std::function<double(double)> fct);