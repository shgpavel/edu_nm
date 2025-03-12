/* SPDX-License-Identifier: Apache-2.0 */

/* Matrix multiplication implementation */

#include <immintrin.h>
#include <openblas/cblas.h>
#include <pthread.h>
#include <unistd.h>

#include "matrix.h"
#include "mm_ob.h"
#include "vector_avx.h"

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

/*
 *
 * function assumes
 * b is transposed,
 * size_rows(a) = size_rows(b) = size_rows(c),
 * size_cols(a) = size_cols(b) = size_cols(c)
 *
 * TODO
 * any matrix
 *
 */

void matrix_mult_fast(matrix *a, matrix *b, matrix *c) {
#ifdef __AVX512CD__
	size_t const block = 4;
#else
	size_t const block = 8;
#endif

	size_t const binr = a->rows / block;
	size_t const bblock = block * block;
	size_t const rowsft = binr * block;
	/* extra row shift */
	size_t inksft = 0;

	if (a->cols != b->rows || a->rows != b->cols || c->rows != a->rows ||
	    c->cols != b->cols || a->rows % block != 0 || a->cols % block != 0) {
		return;
	}

	for (size_t k = 0; k < binr * binr; ++k) {
		if (k % binr == 0 && k != 0) {
			inksft += block * (block - 1) * binr;
		}
		for (size_t h = 0; h < binr; ++h) {
			/* 4x4 block mult */
			for (size_t i = 0; i < block; ++i) {
				for (size_t j = 0; j < block; ++j) {
					size_t left_index =
					    rowsft * i + h * block;
					size_t right_index =
					    rowsft * j + h * block;

					if (k % binr == 0 && k != 0) {
						left_index += bblock * binr;
					} else if (k != 0) {
						right_index += bblock * binr;
					}

#ifdef __AVX512CD__
					__m512d left_row = _mm512_loadu_pd(
					    &a->data->data[left_index]);
					__m512d right_row = _mm512_loadu_pd(
					    &b->data->data[right_index]);

					__m512d mult =
					    _mm512_mul_pd(left_row, right_row);

					matrix_direct(c, i * c->cols + j +
												k * block +
												inksft) += avxreg_sum512(mult);
#else
					__m256d left_row = _mm256_load_pd(
					    &a->data->data[left_index]);
					__m256d right_row = _mm256_load_pd(
					    &b->data->data[right_index]);

					__m256d mult =
					    _mm256_mul_pd(left_row, right_row);

					matrix_direct(
					    c, i * c->cols + j + k * block +
					           inksft) += avxreg_sum(mult);
#endif
				}
			}
		}
	}
}
