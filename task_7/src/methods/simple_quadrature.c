/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Contains quadrature formulas
 * - midpoint (middle rectangles)
 * - left rectangles 
 * - trapezoidal
 * - Simpson's
 */

#include <stdio.h>


#include "../funcs/func.h"
#include "../include/vector.h"

double quad_left(vector *points) {

	double res = 0.0;
	double h = vector_val(points, 1) - vector_val(points, 0);
	for (size_t i = 0; i < points->size - 1; ++i) {
		res += func(vector_val(points, i));
	}
	
	return h * res;
}

double quad_mid(vector *points) {

	double res = 0.0;
	double h = vector_val(points, 1) - vector_val(points, 0);
	for (size_t i = 0; i < points->size - 1; ++i) {
		res += func(vector_val(points, i+1) - h/2);
	}
	
	return h * res;
}

double quad_trapz(vector *points) {

	double res = 0.0;
	double h = vector_val(points, 1) - vector_val(points, 0);
	for (size_t i = 0; i < points->size - 1; ++i) {
		res += func(vector_val(points, i)) + func(vector_val(points, i+1));
	}
	
	return res * h/2;
}

double quad_simpson(vector *points) {

	double res = 0.0;
	double h = vector_val(points, 1) - vector_val(points, 0);
	for (size_t i = 0; i < points->size - 2; i += 2) {
		res += func(vector_val(points, i))
			+ 4 * func(vector_val(points, i+1))
			+ func(vector_val(points, i+2));
	}
	
	return res * h/3;
}
