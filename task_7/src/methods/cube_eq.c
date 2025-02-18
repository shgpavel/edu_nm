/* SPDX-License-Identifier: Apache-2.0 */

#include <jemalloc/jemalloc.h>
#include <math.h>
#include <stdio.h>

#include "polynoms.h"
#include "vector.h"

double pow_s(double x, size_t i) {
	if (i == 2)
		return x * x;
	else if (i == 3)
		return x * x * x;
	else
		return 0.0;
}

vector cube_eq_roots(vector *poly) {
	double a, b, c, d;
	d = vector_val(poly, 0);
	c = vector_val(poly, 1);
	b = vector_val(poly, 2);
	a = vector_val(poly, 3);

	/* y, 2q, 2py, y^3 */
	vector new_var_parts[4];
	for (size_t i = 0; i < 4; ++i) {
		vector_init(&new_var_parts[i], 2, sizeof(double));
		vector_fill_smth(&new_var_parts[i], 0.0);
	}

	vector_val(&new_var_parts[0], 0) = b / 3 * a;
	vector_val(&new_var_parts[0], 1) = 1;

	vector_val(&new_var_parts[1], 0) = 2 * pow_s(b, 3) / 27 * pow_s(a, 3) -
	                                   (b * c / 3 * pow_s(a, 2)) + d / a;

	vector_val(&new_var_parts[2], 0) =
	    2 * ((3 * a * c - pow_s(b, 2)) / 9 * pow_s(a, 2));

	double p = vector_val(&new_var_parts[2], 0) / 2;
	double q = vector_val(&new_var_parts[1], 0) / 2;
	double y_shift = vector_val(&new_var_parts[0], 0);

	vector twopeie = poly_mult(2, &new_var_parts[2], &new_var_parts[0]);
	vector_swap_eff(&new_var_parts[2], &twopeie);

	vector iecubed = poly_mult(3, &new_var_parts[0], &new_var_parts[0],
	                           &new_var_parts[0]);
	vector_swap_eff(&new_var_parts[3], &iecubed);

	vector_free(&twopeie);
	vector_free(&iecubed);

	vector var_sub = poly_sum(3, &new_var_parts[1], &new_var_parts[2],
	                          &new_var_parts[3]);

	for (size_t i = 0; i < 4; ++i) {
		vector_free(&new_var_parts[i]);
	}

	if (pow_s(q, 2) + pow_s(p, 3) >= 0 || p > 0) {
		fprintf(stderr, "Error: it's 3 roots only if q^2 + p^3 < 0\n");

		vector_free(&var_sub);

		vector iamnull;
		vector_init(&iamnull, 0, 0);
		return iamnull;
	}

	double r = q < 0 ? -sqrt(fabs(p)) : sqrt(fabs(p));
	double phi = acos(q / pow_s(r, 3));

	vector roots;
	vector_init(&roots, 3, sizeof(double));
	vector_fill_smth(&roots, 0.0);

	vector_val(&roots, 0) = -2 * r * cos(phi / 3);
	vector_val(&roots, 1) = 2 * r * cos(M_PI / 3 - phi / 3);
	vector_val(&roots, 2) = 2 * r * cos(M_PI / 3 + phi / 3);

	for (size_t i = 0; i < 3; ++i) {
		vector_val(&roots, i) -= y_shift;
	}

	vector_free(&var_sub);

	return roots;
}
