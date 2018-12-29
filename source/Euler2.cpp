#include <cinttypes>

#include "Parameters.hpp"
#include "Derivative.hpp"
#include "Euler1.hpp"

State euler2(const State &r, const Real &h, const Parameters &par)
{
	State sol;
	State s = derivative(euler1(r,0.5*h,par),par);
	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+h*s[i]);
	return sol;
}
