/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Contains advanced quadrature formulas with 3-point interpolation.
 * - Newton-Cotes
 * - Gauss
 */

#include <stdio.h>
#include <tgmath.h>

#include "func.h"
#include "vector.h"
#include "cube_eq.h"
#include "qr_avx.h"

/* Wrong formulas XD */
double calc_moment(size_t k) {
	double inter_len = B - A;
	double moment = pow(inter_len, (double)k + 1 - ALPHA - BETA)
		* tgamma(k + 1 - ALPHA) * tgamma(1 - BETA)
		/ tgamma(k + 2 - ALPHA - BETA);
	return moment;
}

double newton_cotes(vector *points, size_t index) {

	if (points == NULL || index > points->size - 2) {
		return 0.0;
	}
	
	double res = 0.0;

	double moments[3];
	for (size_t i = 0; i < 3; ++i) {
		moments[i] = calc_moment(i);
	}

	double z1 = vector_val(points, index);
	double z12 = (vector_val(points, index + 1) + vector_val(points, index)) / 2;
	double z2 = vector_val(points, index + 1);

	double a1 = moments[2] - moments[1] * (z12 + z2) + moments[0] * z12 * z2;
	a1 /= (z12 - z1) * (z2 - z1);

	double a2 = moments[2] - moments[1] * (z1 + z2) + moments[0] * z1 * z2;
	a2 /= (z12 - z1) * (z2 - z12);
	a2 = -a2;

	double a3 = moments[2] - moments[1] * (z12 + z1) + moments[0] * z12 * z1;
	a3 /= (z2 - z12) * (z2 - z1);

	res = a1 * func(z1) + a2 * func(z12) + a3 * func(z2);

	return res;
}
