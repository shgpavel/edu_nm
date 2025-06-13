/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Contains quadrature formulas
 * - midpoint (middle rectangles)
 * - left rectangles
 * - trapezoidal
 * - Simpson's
 */

#include <stdio.h>

#include "func.h"

double quad_left(double a, double b) {
	return (b - a) * func(a);
}

double quad_mid(double a, double b) {
	return (b - a) * func((a + b) / 2);
}

double quad_trapz(double a, double b) {
	return (func(a) + func(b)) * (b - a) / 2;
}

double quad_simpson(double a, double b) {
	return (func(a) + 4 * func((a + b) / 2) + func(b)) * (b - a) / 6;
}
