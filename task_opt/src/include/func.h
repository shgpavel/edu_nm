/* SPDX-License-Identifier: Apache-2.0 */

#ifndef FUNC_H
#define FUNC_H

#include <math.h>

static inline double func(double x) {
	return exp(-2.6 * x) + exp(0.6 * x) + 1.7 * x;
}

#endif
