#ifndef RK4_H
#define RK4_H

#include <vector>
#include "parametry.hpp"
#include <iostream>
#include <cinttypes>
#include <cmath>
#include "rhs.hpp"

State rk4(const State &r, const Real &h, const Parametry &par);

#endif