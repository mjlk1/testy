#ifndef RK4_H
#define RK4_H

#include <vector>
#include "Parametry.hpp"
#include <iostream>
#include <cinttypes>
#include <cmath>
#include "Derivative.hpp"

State rk4(const State &r, const Real &h, const Parametry &par);

#endif
