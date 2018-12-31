#include <cinttypes>

#include "Parameters.hpp"
#include "Derivative.hpp"

State euler1(const State &r, const Real &h, const Parameters &par)
{
	State sol;
	State s = derivative(r,par);
	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+h*s[i]);
	return sol;
}

