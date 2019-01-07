#include <cmath>

#include "Parameters.hpp"

Real energy(const State &r, const Parameters &par)
{
	Real E;
	E = 0.5*(r[2]*r[2]+r[3]*r[3])-par.g/(sqrt((r[0]-1.0)*(r[0]-1.0)+(r[1]-1.0)*(r[1]-1.0)));
	return E;
}
