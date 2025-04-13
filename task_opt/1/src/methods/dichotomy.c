/* SPDX-License-Identifier: Apache-2.0 */

#include <math.h>
#include <stdio.h>

#include "func.h"

double dichotomy_method(double a, double b, double epsilon) {
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}

	while (b - a > epsilon) {
		double mid = (a + b) * 0.5;
		double f1 = func(mid - epsilon);
		double f2 = func(mid + epsilon);
		if (f1 < f2)
			b = mid;
		else
			a = mid;
	}
	
	return (a + b) * 0.5;
}
