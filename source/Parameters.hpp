#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <cinttypes>

typedef float Real;

typedef std::vector<Real> State; //wektor stanu (x,y,vx,vy)

struct Parameters
{
	Real g, C, n;
};

const int_fast32_t  dist = 2*sizeof(Real)+3;

#endif
