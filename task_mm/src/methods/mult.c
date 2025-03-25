/* SPDX-License-Identifier: Apache-2.0 */

// Matrix multiplication implementations

#include "mult.h"

#include <immintrin.h>
#include <mkl.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "matrix.h"
#include "vector_avx.h"

// B is transposed for all funcs

void matrix_mult_mkl(matrix const* a, matrix const* b, matrix const* c) {
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a->rows, a->rows,
	            a->rows, 1, a->data, a->rows, b->data, a->rows, 0, c->data,
	            a->rows);
}

void matrix_mult_naive(matrix const* a, matrix const* b, matrix const* c) {
	if (a->cols != b->rows || c->rows != a->rows || c->cols != b->cols) {
		return;
	}

	for (size_t i = 0; i < a->rows; ++i) {
		for (size_t j = 0; j < b->cols; ++j) {
			double sum = 0.0;
			for (size_t k = 0; k < a->cols; ++k) {
				sum += matrix_val(a, i, k) * matrix_val(b, j, k);
			}
			matrix_val(c, i, j) = sum;
		}
	}
}

/*
 * size_rows(a) = size_rows(b) = size_rows(c),
 * size_cols(a) = size_cols(b) = size_cols(c)
 * and multiple of 4 or 8
 */

void matrix_mult_fast(matrix const* a, matrix const* b, matrix const* c) {
	size_t const block = 4;

	// beta = 0
	memset(c->data, 0, c->rows * c->cols * sizeof(double));

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
			// multiply a block
			size_t li = h * block + lftsft;
			size_t ri = h * block + risft;

			__m256d cr[4] = {
			    _mm256_load_pd(&c->data[k * block + csft]),
			    _mm256_load_pd(&c->data[c->cols + k * block + csft]),
			    _mm256_load_pd(&c->data[2 * c->cols + k * block + csft]),
			    _mm256_load_pd(
			        &c->data[3 * c->cols + k * block + csft])};

			__m256d left_row[4] = {
			    _mm256_load_pd(&a->data[li]),
			    _mm256_load_pd(&a->data[li + rowsft]),
			    _mm256_load_pd(&a->data[li + 2 * rowsft]),
			    _mm256_load_pd(&a->data[li + 3 * rowsft])};

			__m256d right_row[4] = {
			    _mm256_load_pd(&b->data[ri]),
			    _mm256_load_pd(&b->data[rowsft + ri]),
			    _mm256_load_pd(&b->data[2 * rowsft + ri]),
			    _mm256_load_pd(&b->data[3 * rowsft + ri])};

			__m256d add[4] = {
			    _mm256_set_pd(
			        avxreg_sum(_mm256_mul_pd(left_row[0], right_row[3])),
			        avxreg_sum(_mm256_mul_pd(left_row[0], right_row[2])),
			        avxreg_sum(_mm256_mul_pd(left_row[0], right_row[1])),
			        avxreg_sum(
			            _mm256_mul_pd(left_row[0], right_row[0]))),
			    _mm256_set_pd(
			        avxreg_sum(_mm256_mul_pd(left_row[1], right_row[3])),
			        avxreg_sum(_mm256_mul_pd(left_row[1], right_row[2])),
			        avxreg_sum(_mm256_mul_pd(left_row[1], right_row[1])),
			        avxreg_sum(
			            _mm256_mul_pd(left_row[1], right_row[0]))),
			    _mm256_set_pd(
			        avxreg_sum(_mm256_mul_pd(left_row[2], right_row[3])),
			        avxreg_sum(_mm256_mul_pd(left_row[2], right_row[2])),
			        avxreg_sum(_mm256_mul_pd(left_row[2], right_row[1])),
			        avxreg_sum(
			            _mm256_mul_pd(left_row[2], right_row[0]))),
			    _mm256_set_pd(
			        avxreg_sum(_mm256_mul_pd(left_row[3], right_row[3])),
			        avxreg_sum(_mm256_mul_pd(left_row[3], right_row[2])),
			        avxreg_sum(_mm256_mul_pd(left_row[3], right_row[1])),
			        avxreg_sum(
			            _mm256_mul_pd(left_row[3], right_row[0])))};

			_mm256_store_pd(&c->data[k * block + csft],
			                _mm256_add_pd(add[0], cr[0]));
			_mm256_store_pd(&c->data[c->cols + k * block + csft],
			                _mm256_add_pd(add[1], cr[1]));
			_mm256_store_pd(&c->data[2 * c->cols + k * block + csft],
			                _mm256_add_pd(add[2], cr[2]));
			_mm256_store_pd(&c->data[3 * c->cols + k * block + csft],
			                _mm256_add_pd(add[3], cr[3]));
		}
	}
}
