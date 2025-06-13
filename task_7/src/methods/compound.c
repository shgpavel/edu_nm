/* SPDX-License-Identifier: Apache-2.0 */

#include "compound.h"

#include <stdio.h>
#include <stdlib.h>

#include "advanced_quadrature.h"
#include "func.h"

double compound(size_t n, quad_func qf) {
	double step = (B - A) / n;
	double total = 0.0;

	for (size_t i = 0; i < n; ++i) {
		double ai = A + i * step;
		double bi = ai + step;
		total += qf(ai, bi);
	}
	return total;
}

/* Mixed some-point Newton-Cotes and 3-p Gauss */
double compound_ncwgau(size_t n, size_t nc_count, quad_func_any_n nc,
                       quad_func gs) {
	double step = (B - A) / n;
	double total = 0.0;

	for (size_t i = 0; i < n; ++i) {
		double ai = A + i * step;
		double bi = ai + step;

		if (rand() % 2 == 0) {
			total += nc(nc_count, ai, bi);
		} else {
			total += gs(ai, bi);
		}
	}
	return total;
}
