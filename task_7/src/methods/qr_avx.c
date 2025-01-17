/* SPDX-License-Identifier: Apache-2.0 */

/* It is not finished */
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "matrix.h"
#include "vector.h"
#include "vector_avx.h"

vector *qr_avx(matrix *a, vector *b) {
	vector *sol = (vector *)malloc(sizeof(vector));

	matrix q;
	matrix r;
	vector y;
	vector ys;

	matrix_cinit(&q, a->rows, a->cols, sizeof(double));
	matrix_cinit(&r, a->rows, a->cols, sizeof(double));
	vector_cinit(&y, b->size, sizeof(double));
	vector_cinit(&ys, b->size, sizeof(double));
	vector_cinit(sol, b->size, sizeof(double));

	for (size_t i = 0; i < a->cols - 1; ++i) {
		for (size_t j = i; j < a->rows; ++j) {
			vector_val(&y, j) = matrix_val(a, j, i);
		}

		vector_val(&y, i) = vector_val(&y, i) - vector_norm2_avx(&y, i);

		double y_norm2 = vector_norm2_avx(&y, i);
		for (size_t j = i; j < a->rows; ++j) {
			vector_val(&y, j) = vector_val(&y, j) / y_norm2;
		}

		for (size_t j = i; j < a->rows; ++j) {
			for (size_t k = i; k < a->cols; ++k) {
				matrix_val(&q, j, k) = -2.0 *
				                       vector_val(&y, j) *
				                       vector_val(&y, k);
			}
			matrix_val(&q, j, j) = matrix_val(&q, j, j) + 1;
		}

		for (size_t j = i; j < a->rows; ++j) {
			for (size_t k = i; k < a->cols; ++k) {
				double sum = 0.0;
				for (size_t l = i; l < a->rows; ++l) {
					sum = sum + (matrix_val(&q, j, l) *
					             matrix_val(a, l, k));
				}
				matrix_val(&r, j, k) = sum;
			}
		}

		for (size_t j = i; j < a->rows; ++j) {
			for (size_t k = i; k < a->cols; ++k) {
				matrix_val(a, j, k) = matrix_val(&r, j, k);
			}
		}

		for (size_t j = i; j < a->rows; ++j) {
			double sum = 0.0;
			for (size_t k = i; k < a->cols; ++k) {
				sum = sum +
				      matrix_val(&q, j, k) * vector_val(b, k);
			}
			vector_val(&ys, j) = sum;
		}

		for (size_t j = i; j < a->rows; ++j) {
			vector_val(b, j) = vector_val(&ys, j);
		}
	}

	for (size_t i = a->rows; i > 0;) {
		--i;
		vector_val(sol, i) = vector_val(b, i);
		for (size_t j = a->rows; j > i + 1;) {
			--j;
			vector_val(sol, i) =
			    vector_val(sol, i) -
			    matrix_val(a, i, j) * vector_val(sol, j);
		}
		vector_val(sol, i) = vector_val(sol, i) / matrix_val(a, i, i);
	}

	matrix_free(&q);
	matrix_free(&r);
	vector_free(&y);
	vector_free(&ys);

	return sol;
}
