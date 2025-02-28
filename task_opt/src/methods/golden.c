/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>
#include <math.h>

#include "func.h"

#define DEBUG_INFO

double golden(double a, double b, double epsilon) {
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}

	size_t itr = 0;
	
	const double phi = (1 + sqrt(5)) / 2;
	const double inv = 1 / phi;

	double x1 = b - inv * (b - a);
	double x2 = a + inv * (b - a);

	double f1 = func(x1);
	double f2 = func(x2);

	while ((b - a) > epsilon) {
		if (f1 < f2) {
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = b - inv * (b - a);
			f1 = func(x1);
		} else {
			a = x1;
			x1 = x2;
			f1 = f2;
			x2 = a + inv * (b - a);
			f2 = func(x2);
		}
	}

	return (a + b) / 2;
}
