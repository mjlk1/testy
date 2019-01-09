#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <cinttypes>
#include <immintrin.h>
#include <stdio.h>

typedef double Real;

typedef std::vector<Real> State;

typedef __m256d StateAVX; 

struct Parameters
{
	Real g, C, n;
};

const int_fast32_t  dist = 2*sizeof(Real)+3;

#endif
