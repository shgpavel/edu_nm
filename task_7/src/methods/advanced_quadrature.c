/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Contains advanced quadrature formulas with 3-point interpolation.
 * - Newton-Cotes
 * - Gauss
 */

#include <math.h>
#include <stdio.h>

#include "../funcs/func.h"
#include "../include/vector.h"
#include "../include/cube_eq.h"

double newton_cotes(vector *points, size_t index) {

	if (points == NULL || index > points->size - 2) {
		return 0.0;
	}
	
	double res = 0.0;
	double z1 = vector_val(points, index);
	double z12 = (vector_val(points, index + 1) + vector_val(points, index)) / 2;
	double z2 = vector_val(points, index + 1);

	vector moments;
	vector_init(&moments, 3, sizeof(double));
	vector_fill_smth(&moments, 0.0);
	
	vector_val(&moments, 0) = (pow(z2-A, 1-ALPHA) - pow(z1-A, 1-ALPHA)) / (1 - ALPHA);
	vector_val(&moments, 1) = (pow(z2-A, 2-ALPHA) - pow(z1-A, 2-ALPHA)) / (2 - ALPHA)
		+ A * vector_val(&moments, 0);
	vector_val(&moments, 2) = (pow(z2-A, 3-ALPHA) - pow(z1-A, 3-ALPHA)) / (3 - ALPHA)
		+ 2 * A * vector_val(&moments, 1) - A*A * vector_val(&moments, 0);

	/*
	vector xminus;
	vector_init(&xminus, 2, sizeof(double));
	vector_val(&xminus, 0) = -A;
	vector_val(&xminus, 1) = 1;
	
	for (size_t i = 0; i < points->size - 1; ++i) {
		res += func(vector_val(points, i));
	}
	*/
	
	return res;
}
