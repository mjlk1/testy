#ifndef PARAMETRY_H
#define PARAMETRY_H

using namespace std;
#include <iostream>
#include <cinttypes>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstring>

typedef float Real;
typedef vector<Real> State; //wektor stanu (x,y,vx,vy)
struct Parametry
{
	Real g,C;
};

const Real mass = 1.0;
const Real n = 2.0;
const int_fast32_t  dist = 2*sizeof(Real)+3;

#endif
