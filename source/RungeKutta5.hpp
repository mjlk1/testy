#ifndef RUNGE_KUTTA_5_HPP
#define RUNGE_KUTTA_5_HPP

#include "Parameters.hpp"

State rk5(const State &r, const Real &h, const Parameters &par);

#endif