#include <vector>
#include "parametry.hpp"
#include <iostream>
#include <cinttypes>
#include <cmath>
#include "rhs.hpp"

State rhs(const State &r, const Parametry &par)
{
	State deriv;
	Real l = -par.C*pow(r[2]*r[2]+r[3]*r[3], (n-1)/2)/mass;
	deriv.push_back(r[2]);
	deriv.push_back(r[3]);
	deriv.push_back(l*r[2]);
	deriv.push_back(l*r[3]+par.g);
	return deriv;
}
