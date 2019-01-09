#ifndef RUNGE_KUTTA_5_AVX_HPP
#define RUNGE_KUTTA_5_AVX_HPP

#include "Parameters.hpp"
#include "DerivativeAVX.hpp"

StateAVX rk5AVX(const StateAVX &r, const Real &h, const Parameters &par);

inline StateAVX rk5AVXInline(const StateAVX &r, const Real &h, const Parameters &par)
{
	StateAVX k1, k2, k3, k4, k5, k6, s1, s2, s3, s4, s5, s6;
	StateAVX H = _mm256_set1_pd(h);

	StateAVX b[7][6];
	StateAVX c[7];
	StateAVX c2[7];

	// b i c to wspolczynniki z RK5
	b[2][1] = _mm256_set1_pd(0.2);

	b[3][1] = _mm256_set1_pd(3.0/40.0);
	b[3][2] = _mm256_set1_pd(9.0/40.0);

	b[4][1] = _mm256_set1_pd(0.3);
	b[4][2] = _mm256_set1_pd(-0.9);
	b[4][3] = _mm256_set1_pd(6.0/5.0);

	b[5][1] = _mm256_set1_pd(-11.0/54.0);
	b[5][2] = _mm256_set1_pd(2.5);
	b[5][3] = _mm256_set1_pd(-70.0/27.0);
	b[5][4] = _mm256_set1_pd(35.0/27.0);

	b[6][1] = _mm256_set1_pd(1631.0/55296.0);
	b[6][2] = _mm256_set1_pd(175.0/512.0);
	b[6][3] = _mm256_set1_pd(575.0/13824.0);
	b[6][4] = _mm256_set1_pd(44275.0/110592.0);
	b[6][5] = _mm256_set1_pd(253.0/4096.0);

	c[1] = _mm256_set1_pd(37.0/378.0);
	c[2] = _mm256_set1_pd(0.0);
	c[3] = _mm256_set1_pd(250.0/621.0);
	c[4] = _mm256_set1_pd(125.0/594.0);
	c[5] = _mm256_set1_pd(0.0);
	c[6] = _mm256_set1_pd(512.0/1771.0);

	s1 = derivativeAVXInline(r,par);
	k1 = H*s1;

	s2 = derivativeAVXInline(r+b[2][1]*k1,par);
	k2 = H*s2;

	s3 = derivativeAVXInline(r+b[3][1]*k1+b[3][2]*k2,par);
	k3 = H*s3;

	s4 = derivativeAVXInline(r+b[4][1]*k1+b[4][2]*k2+b[4][3]*k3,par);
	k4 = H*s4;

	s5 = derivativeAVXInline(r+b[5][1]*k1+b[5][2]*k2+b[5][3]*k3+b[5][4]*k4,par);
	k5 = H*s5;

	s6 = derivativeAVXInline(r+b[6][1]*k1+b[6][2]*k2+b[6][3]*k3+b[6][4]*k4+b[6][5]*k5,par);
	k6 = H*s6;

	return r+c[1]*k1+c[3]*k3+c[4]*k4+c[6]*k6;
}

#endif
