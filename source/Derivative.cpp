#include <cmath>

#include "Parameters.hpp"

State derivative(const State &r, const Parameters &par)
{
	State deriv;
	Real k = -par.g*(par.C+1.0)/pow((r[0])*(r[0])+(r[1])*(r[1]),(par.C+3.0)/2.0);
	deriv.push_back(r[2]);
	deriv.push_back(r[3]);
	deriv.push_back(k*(r[0]));
	deriv.push_back(k*(r[1]));
	return deriv;
}
