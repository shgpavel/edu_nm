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

#ifdef DEBUG_INFO
	size_t itr = 1;
#endif

	double x3, f3;

	while (1) {
		x3 = x2 + step;
		f3 = func(x3);

#ifdef DEBUG_INFO
		if (itr == 1)
			printf(
			    "step 1\n"
			    "i |     x_2     |     x_3     |    f(x_2)   |"
			    "     f(x_3)  |  step\n");

		char pool[5][25];
		double values[5] = {x2, x3, f2, f3, step};

		for (size_t i = 0; i < 5; ++i)
			format_hex(pool[i], 25, values[i]);

		printf("%zu | %s | %s | %s | %s | %s\n", itr, pool[0], pool[1],
		       pool[2], pool[3], pool[4]);

		++itr;
#endif

		if (f3 > f2) break;

		x1 = x2;
		x2 = x3;
		f1 = f2;
		f2 = f3;
	}

#ifdef DEBUG_INFO
	itr = 1;
#endif

	int const n = -5;
	double const epsilon = 8 * pow(16, -4) * pow(2, n);

	int change = 1;
	while (fabs(step) >= epsilon) {
		if (change) step /= 2;

		double x_temp = x2 + step;
		double f_temp = func(x_temp);

		if (f2 > f_temp) {
			x1 = x2;
			f1 = f2;
			x2 = x_temp;
			f2 = f_temp;
			change = 1;
		} else {
			if (!change) {
				x3 = x_temp;
				f3 = f_temp;
				change = 1;
			} else {
				x3 = x1;
				f3 = f1;
				x1 = x_temp;
				f1 = f_temp;
				step = -step;
				change = 0;
			}
		}

#ifdef DEBUG_INFO
		if (itr == 1)
			printf(
			    "\nstep 2\n"
			    "i |     x_1      |     x_2     |     x_3     |"
			    "    f(x_1)   |"
			    "    f(x_2)   |   f(x_3)   |   step\n");

		char pool[7][25];
		double values[7] = {x1, x2, x3, f1, f2, f3, step};

		for (size_t i = 0; i < 7; ++i)
			format_hex(pool[i], 25, values[i]);

		printf("%zu | %s | %s | %s | %s | %s | %s | %s\n", itr, pool[0],
		       pool[1], pool[2], pool[3], pool[4], pool[5], pool[6]);
		++itr;
#endif
	}

	return x2;
}
