#include <vector>
#include "parametry.hpp"
#include <iostream>
#include <cinttypes>
#include <cmath>
#include "rhs.hpp"

State rk4(const State &r, const Real &h, const Parametry &par)
{
	State k1, k2, k3, k4, yk1, yk2, yk3, sol;

	State s = rhs(r,par);
	for (int_fast32_t i=0;i<4;++i)
		k1.push_back(h*s[i]); //obliczamy k1

	for (int_fast32_t i=0;i<4;++i)
		yk1.push_back(r[i]+0.5*k1[i]); //obliczamy yk1

	State sk1 = rhs(yk1,par);
	for (int_fast32_t i=0;i<4;++i)
		k2.push_back(h*sk1[i]); //obliczamy k2

	for (int_fast32_t i=0;i<4;++i)
		yk2.push_back(r[i]+0.5*k2[i]); //obliczamy yk2

	State sk2 = rhs(yk2,par);
	for (int_fast32_t i=0;i<4;++i)
		k3.push_back(h*sk2[i]); //obliczamy k3

	for (int_fast32_t i=0;i<4;++i)
		yk3.push_back(r[i]+k3[i]); //obliczamy yk3

	State sk3 = rhs(yk3,par);
	for (int_fast32_t i=0;i<4;++i)
		k4.push_back(h*sk3[i]); //obliczamy k3

	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6); //obliczamy k3

	return sol;
}
