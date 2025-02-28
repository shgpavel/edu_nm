/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>
#include <math.h>

#include "func.h"

#define DEBUG_INFO

double dichotomy_method(double a, double b, double epsilon, double delta) {
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}

	delta = fmin(delta, (b - a) * 0.49);

	while ((b - a) > epsilon) {
		double mid = (a + b) / 2;
		double x1 = mid - delta / 2;
		double x2 = mid + delta / 2;

		double f1 = func(x1);
		double f2 = func(x2);

		if (f1 < f2) {
			b = x2;
		} else {
			a = x1;
		}
		
	}

	return (a + b) / 2;
}
