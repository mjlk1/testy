#include <cinttypes>

#include "Parameters.hpp"
#include "Derivative.hpp"

State rk5(const State &r, const Real &h, const Parameters &par)
{
	State k1, k2, k3, k4, k5, k6, yk1, yk2, yk3, yk4, yk5, sol;

	State s = derivative(r,par);
	for (int_fast32_t i=0;i<4;++i)
		k1.push_back(h*s[i]); //obliczamy k1

	for (int_fast32_t i=0;i<4;++i)
		yk1.push_back(r[i]+0.2*k1[i]); //obliczamy yk1

	State sk1 = derivative(yk1,par);
	for (int_fast32_t i=0;i<4;++i)
		k2.push_back(h*sk1[i]); //obliczamy k2

	for (int_fast32_t i=0;i<4;++i)
		yk2.push_back(r[i]+(3.0/40.0)*k1[i]+(9.0/40.0)*k2[i]); //obliczamy yk2

	State sk2 = derivative(yk2,par);
	for (int_fast32_t i=0;i<4;++i)
		k3.push_back(h*sk2[i]); //obliczamy k3

	for (int_fast32_t i=0;i<4;++i)
		yk3.push_back(r[i]+0.3*k1[i]-0.9*k2[i]+(6.0/5.0)*k3[i]); //obliczamy yk3

	State sk3 = derivative(yk3,par);
	for (int_fast32_t i=0;i<4;++i)
		k4.push_back(h*sk3[i]); //obliczamy k4

	for (int_fast32_t i=0;i<4;++i)
		yk4.push_back(r[i]-(11.0/54.0)*k1[i]+2.5*k2[i]-(70.0/27.0)*k3[i]+(35.0/27.0)*k4[i]); //obliczamy yk4

	State sk4 = derivative(yk4,par);
	for (int_fast32_t i=0;i<4;++i)
		k5.push_back(h*sk4[i]); //obliczamy k5

	for (int_fast32_t i=0;i<4;++i)
		yk5.push_back(r[i]+(1631.0/55296.0)*k1[i]+(175.0/512.0)*k2[i]+(575.0/13824.0)*k3[i]+(44275.0/110592.0)*k4[i]+(253.0/4096.0)*k5[i]); //obliczamy yk5

	State sk5 = derivative(yk5,par);
	for (int_fast32_t i=0;i<4;++i)
		k6.push_back(h*sk5[i]); //obliczamy k6

	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+(37.0/378.0)*k1[i]+(250.0/621.0)*k3[i]+(125.0/594.0)*k4[i]+(512.0/1771.0)*k6[i]); 

	return sol;
}
