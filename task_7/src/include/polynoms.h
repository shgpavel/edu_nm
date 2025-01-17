/* SPDX-License-Identifier: Apache-2.0 */

#ifndef POLYNOMS_H
#define POLYNOMS_H

#include "vector.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

vector poly_mult(size_t, ...);
vector poly_sum(size_t, ...);
double poly_val(vector *, double);

#endif
