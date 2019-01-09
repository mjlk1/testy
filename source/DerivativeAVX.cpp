#include <cmath>
#include <immintrin.h>

#include "Parameters.hpp"

StateAVX derivativeAVX(const StateAVX &r, const Parameters &par)
{
	Real k = -par.g*(par.C+1.0)/pow(r[0]*r[0]+r[1]*r[1],(par.C+3.0)/2.0);
	StateAVX deriv;
	deriv[0] = r[2];
	deriv[1] = r[3];
	deriv[2] = k*r[0];
	deriv[3] = k*r[1];

	return deriv;
}
