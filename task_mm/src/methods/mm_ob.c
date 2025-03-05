/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Matrix multiplication
 */

#include <emmintrin.h>
#include <immintrin.h>
#include <openblas/cblas.h>
#include <vector_avx.h>

#include "matrix.h"

void matrix_mult_openblas(matrix *a, matrix *b, matrix *c) {
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a->rows, a->rows,
	            a->rows, 1, a->data->data, a->rows, b->data->data, a->rows,
	            0, c->data->data, a->rows);
}

void matrix_mult_naive(matrix *a, matrix *b, matrix *c) {
	if (a->cols != b->rows || c->rows != a->rows || c->cols != b->cols) {
		return;
	}

	for (size_t i = 0; i < a->rows; ++i) {
		for (size_t j = 0; j < b->cols; ++j) {
			double sum = 0.0;
			for (size_t k = 0; k < a->cols; ++k) {
				sum +=
				    matrix_val(a, i, k) * matrix_val(b, k, j);
			}
			matrix_val(c, i, j) = sum;
		}
	}
}

/* TODO
 * any matrix
 * function assumes b is transposed
 */
void matrix_mult_fast(matrix *a, matrix *b, matrix *c) {
	if (a->cols != b->rows || a->rows != b->cols || c->rows != a->rows ||
	    c->cols != b->cols || a->rows % 4 != 0 || a->cols % 4 != 0) {
		return;
	}

	for (size_t i = 0; i < a->rows; ++i) {
		for (size_t j = 0; j < b->cols; ++j) {
			// for (size_t k = 0; k < 4; ++k) {
			__m256 left_row = _mm256_load_pd(&a->data->data[4 * i]);
			__m256 right_row =
			    _mm256_load_pd(&b->data->data[4 * j]);
			__m256 mult = _mm256_mul_pd(left_row, right_row);

			// print_avxreg(mult);
			// printf("\n");

			__m128d low = _mm256_castpd256_pd128(mult);
			__m128d high = _mm256_extractf128_pd(mult, 1);

			/* TODO find better way to get sum of AVX reg */
			__m128d sum_r = _mm_add_pd(low, high);
			matrix_val(c, i, j) =
			    _mm_cvtsd_f64(sum_r) +
			    _mm_cvtsd_f64(_mm_unpackhi_pd(sum_r, sum_r));
		}
	}
}
