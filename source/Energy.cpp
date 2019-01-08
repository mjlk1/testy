#include <cmath>

#include "Parameters.hpp"

Real energy(const State &r, const Parameters &par)
{
	Real E;
	E = 0.5*(r[2]*r[2]+r[3]*r[3])-par.g/pow(r[0]*r[0]+r[1]*r[1],(par.C+1.0)/2.0);
	return E;
}
