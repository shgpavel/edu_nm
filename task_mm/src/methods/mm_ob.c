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

void matrix_mult_mkl(matrix const* restrict a, matrix const* restrict b,
                     matrix const* restrict c) {
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a->rows, a->rows, a->rows, 1,
	            a->data->data, a->rows, b->data->data, a->rows, 0, c->data->data, a->rows);
}

void matrix_mult_naive(matrix const* restrict a, matrix const* restrict b,
                       matrix const* restrict c) {
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
 *
 * TODO
 * any matrix, better cache util
 */

void matrix_mult_fast(matrix const* restrict a, matrix const* restrict b,
                      matrix const* restrict c) {
#ifdef __AVX512CD__
	size_t const block = 8;
#elifdef __AVX2__
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
			/* multiply a block */
			for (size_t i = 0; i < block; ++i) {
				size_t left_index = rowsft * i + h * block + lftsft;
				size_t ri_base = h * block + risft;

#ifdef __AVX512CD__
				__m512d cp =
				    _mm512_load_pd(&c->data->data[i * c->cols + k * block + csft]);

				__m512d left_row = _mm512_load_pd(&a->data->data[left_index]);

				__m512d right_row_0 = _mm512_load_pd(&b->data->data[ri_base]);
				__m512d right_row_1 =
				    _mm512_load_pd(&b->data->data[rowsft + ri_base]);
				__m512d right_row_2 =
				    _mm512_load_pd(&b->data->data[rowsft * 2 + ri_base]);
				__m512d right_row_3 =
				    _mm512_load_pd(&b->data->data[rowsft * 3 + ri_base]);
				__m512d right_row_4 =
				    _mm512_load_pd(&b->data->data[rowsft * 4 + ri_base]);
				__m512d right_row_5 =
				    _mm512_load_pd(&b->data->data[rowsft * 5 + ri_base]);
				__m512d right_row_6 =
				    _mm512_load_pd(&b->data->data[rowsft * 6 + ri_base]);
				__m512d right_row_7 =
				    _mm512_load_pd(&b->data->data[rowsft * 7 + ri_base]);

				__m512d add = _mm512_set_pd(
				    avxreg_sum512(_mm512_mul_pd(left_row, right_row_7)),
				    avxreg_sum512(_mm512_mul_pd(left_row, right_row_6)),
				    avxreg_sum512(_mm512_mul_pd(left_row, right_row_5)),
				    avxreg_sum512(_mm512_mul_pd(left_row, right_row_4)),
				    avxreg_sum512(_mm512_mul_pd(left_row, right_row_3)),
				    avxreg_sum512(_mm512_mul_pd(left_row, right_row_2)),
				    avxreg_sum512(_mm512_mul_pd(left_row, right_row_1)),
				    avxreg_sum512(_mm512_mul_pd(left_row, right_row_0)));
				__m512 res = _mm512_add_pd(add, cp);
				_mm512_store_pd(&c->data->data[i * c->cols + k * block + csft],
				                res);

#elifdef __AVX2__

				__m256d cp =
				    _mm256_load_pd(&c->data->data[i * c->cols + k * block + csft]);

				__m256d left_row = _mm256_load_pd(&a->data->data[left_index]);

				__m256d right_row_0 = _mm256_load_pd(&b->data->data[ri_base]);
				__m256d right_row_1 =
				    _mm256_load_pd(&b->data->data[rowsft + ri_base]);
				__m256d right_row_2 =
				    _mm256_load_pd(&b->data->data[2 * rowsft + ri_base]);
				__m256d right_row_3 =
				    _mm256_load_pd(&b->data->data[3 * rowsft + ri_base]);

				__m256d add = _mm256_set_pd(
				    avxreg_sum(_mm256_mul_pd(left_row, right_row_3)),
				    avxreg_sum(_mm256_mul_pd(left_row, right_row_2)),
				    avxreg_sum(_mm256_mul_pd(left_row, right_row_1)),
				    avxreg_sum(_mm256_mul_pd(left_row, right_row_0)));
				__m256 res = _mm256_add_pd(add, cp);
				_mm256_store_pd(&c->data->data[i * c->cols + k * block + csft],
				                res);
#endif
			}
		}
	}
}
