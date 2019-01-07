#include <cmath>

#include "Parameters.hpp"

State derivative(const State &r, const Parameters &par)
{
	State deriv;
	Real l = -par.C*pow(r[2]*r[2]+r[3]*r[3], (par.n-1.0)/2.0);
	deriv.push_back(r[2]);
	deriv.push_back(r[3]);
	deriv.push_back(l*r[2]-par.g*(r[0]-1.0)/pow((r[0]-1.0)*(r[0]-1.0)+(r[1]-1.0)*(r[1]-1.0),1.5));
	deriv.push_back(l*r[3]-par.g*(r[1]-1.0)/pow((r[0]-1.0)*(r[0]-1.0)+(r[1]-1.0)*(r[1]-1.0),1.5));
	return deriv;
}
