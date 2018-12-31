#ifndef RUNGE_KUTTA_4_HPP
#define RUNGE_KUTTA_4_HPP

#include "Parameters.hpp"

State rk4(const State &r, const Real &h, const Parameters &par);

#endif
