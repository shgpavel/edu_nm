/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Matrix multiplication implementations
 * All functions work with b transposed
 */

#include "mm_ob.h"

#include <immintrin.h>
#include <mkl.h>
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>

#include "matrix.h"
#include "vector_avx.h"

void matrix_mult_mkl(matrix *a, matrix *b, matrix *c) {
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a->rows, a->rows,
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
				    matrix_val(a, i, k) * matrix_val(b, j, k);
			}
			matrix_val(c, i, j) = sum;
		}
	}
}

/*
 * size_rows(a) = size_rows(b) = size_rows(c),
 * size_cols(a) = size_cols(b) = size_cols(c)
 * and multiple of 4 or 8
 *
 * TODO
 * any matrix, FMA, cache usg
 */

void matrix_mult_fast(matrix *a, matrix *b, matrix *c) {
#ifdef __AVX512CD__
	size_t const block = 8;
#else
	size_t const block = 4;
#endif

	size_t const blinrow = a->rows / block;
	size_t const blocksq = block * block;
	size_t const rowsft = blinrow * block;

	size_t csft = 0;
	size_t lftsft = 0;
	size_t risft = 0;

	for (size_t k = 0; k < blinrow * blinrow; ++k) {
		if (k % blinrow == 0 && k != 0) {
			csft += block * (block - 1) * blinrow;
			lftsft += blocksq * blinrow;
			risft = 0;
		} else if (k != 0) {
			risft += blocksq * blinrow;
		}
		for (size_t h = 0; h < blinrow; ++h) {
			/* block mult */
			for (size_t i = 0; i < block; ++i) {
				for (size_t j = 0; j < block; ++j) {
					size_t left_index =
					    rowsft * i + h * block + lftsft;
					size_t right_index =
					    rowsft * j + h * block + risft;

#ifdef __AVX512CD__
					__m512d left_row = _mm512_load_pd(
					    &a->data->data[left_index]);
					__m512d right_row = _mm512_load_pd(
					    &b->data->data[right_index]);

					__m512d mult =
					    _mm512_mul_pd(left_row, right_row);

					matrix_direct(c, i * c->cols + j +
					                     k * block +
					                     inksft) +=
					    avxreg_sum512(mult);
#else
					__m256d left_row = _mm256_load_pd(
					    &a->data->data[left_index]);
					__m256d right_row = _mm256_load_pd(
					    &b->data->data[right_index]);

					__m256d mult =
					    _mm256_mul_pd(left_row, right_row);

#ifdef DEBUG
					printf(
					    "k = %zu, h = %zu, ij %zu %zu, "
					    "liri %zu %zu\n",
					    k, h, i, j, left_index,
					    right_index);
					avxreg_print(left_row);
					avxreg_print(right_row);
					printf("\n");
#endif

					matrix_direct(
					    c, i * c->cols + j + k * block +
					           csft) += avxreg_sum(mult);
#endif
				}
			}
		}
	}
}
