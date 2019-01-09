#ifndef RUNGE_KUTTA_5_AVX_HPP
#define RUNGE_KUTTA_5_AVX_HPP

#include "Parameters.hpp"

StateAVX rk5AVX(const StateAVX &r, const Real &h, const Parameters &par);

#endif