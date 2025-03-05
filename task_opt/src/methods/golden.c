/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>
#include <math.h>

#include "func.h"

double golden(double a, double b, double epsilon) {
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}

	size_t itr = 1;
	
	const double phi = (1 + sqrt(5)) / 2;
	const double inv = 1 / phi;

	double c = b - inv * (b - a);
	double d = a + inv * (b - a);

	double f1 = func(c);
	double f2 = func(d);

#ifdef DEBUG_INFO
	printf("i: %zu c: %.3a d: %.3a f(c): %.3a f(d-c): %.3a\n",
				 itr, c, d, func(c), func(d-c));
#endif
	
	while ((b - a) > epsilon) {
		if (f1 < f2) {
			b = d;
			d = c;
			f2 = f1;
			c = b - inv * (b - a);
			f1 = func(c);
		} else {
			a = c;
			c = d;
			f1 = f2;
			d = a + inv * (b - a);
			f2 = func(d);
		}
		++itr;
		
#ifdef DEBUG_INFO
		printf("i: %zu c: %.3a d: %.3a f(c): %.3a f(d-c): %.3a\n",
					 itr, c, d, func(c), func(d-c));
#endif
	}
	
	return (a + b) / 2;
}
