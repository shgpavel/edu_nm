/* SPDX-License-Identifier: Apache-2.0 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "func.h"

double direct_search(double init, double step) {
	double x1 = init;
	double x2 = x1 + step;
	double f1 = func(x1);
	double f2 = func(x2);

	if (f1 < f2) {
		step = -step;
		double temp = x1;
		x1 = x2;
		x2 = temp;
		temp = f1;
		f1 = f2;
		f2 = temp;
	}

	while (1) {
		double x3 = x2 + step;
		double f3 = func(x3);

		if (f3 <= f2) {
			x1 = x2;
			x2 = x3;
			f1 = f2;
			f2 = f3;
		} else {
			break;
		}
	}

	int change = 1;
	int n = -5;
	double epsilon = 8 * pow(16, -4) * pow(2, n);

	while (fabs(step) >= epsilon) {
		if (change) {
			step /= 2;
		}

		double x_temp = x2 + step;
		double f_temp = func(x_temp);

		if (f2 > f_temp) {
			x1 = x2;
			f1 = f2;
			x2 = x_temp;
			f2 = f_temp;
			change = 1;
		} else {
			if (change) {
				x1 = x_temp;
				f1 = f_temp;
				step = -step;
				change = 0;
			} else {
				change = 1;
			}
		}
	}

	return x2;
}
