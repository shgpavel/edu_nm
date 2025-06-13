/* SPDX-License-Identifier: Apache-2.0 */

#include <math.h>
#include <stdio.h>

#include "func.h"

double dichotomy_method(double a, double b, double epsilon, double delta) {
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}

#ifdef DEBUG_INFO
	size_t itr = 1;
#endif

	while (b - a > epsilon) {
		double mid = (a + b) * 0.5;

		double c = mid - delta;
		double d = mid + delta;

		double f1 = func(c);
		double f2 = func(d);

		if (f1 < f2)
			b = d;
		else
			a = c;

#ifdef DEBUG_INFO
		if (itr == 1)
			printf(
			    "i |     x_a     |     x_b     |     x_c     |"
			    "     x_d     "
			    "|    f(x_c)   | f(x_d)-f(x_c) \n");

		char pool[6][25];
		double values[6] = {a, b, c, d, f1, f2 - f1};

		for (size_t i = 0; i < 6; ++i)
			format_hex(pool[i], 25, values[i]);

		printf("%zu | %s | %s | %s | %s | %s | %s\n", itr, pool[0],
		       pool[1], pool[2], pool[3], pool[4], pool[5]);

		++itr;
#endif
	}

	return (a + b) * 0.5;
}
