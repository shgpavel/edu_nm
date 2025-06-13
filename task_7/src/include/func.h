/* SPDX-License-Identifier: Apache-2.0 */

#ifndef FUNC_H
#define FUNC_H

#include <math.h>

static double const A = 1.1;
static double const B = 2.3;
static double const ALPHA = 0.8;
static double const BETA = 0;

static inline double func(double x) {
	return 3.5 * cos(0.7 * x) * exp(-5 * x / 3) +
	       2.4 * sin(5.5 * x) * exp(-3 * x / 4) + 5;
}

#endif
