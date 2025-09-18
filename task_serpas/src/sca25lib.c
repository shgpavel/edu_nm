/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Copyright (C) 2025 Pavel Shago <pavel@shago.dev>
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <omp.h>

// gauss low(x86) high(e2k)
double *
funcTask1(double **a, double *y, int N)
{
	omp_set_num_threads(3);

	double *sol = (double *)malloc(N * sizeof(double));

	for (int i = 0; i < N; ++i) {
		int lead = i;

		for (int j = i + 1; j < N; ++j) {
			if (fabs(a[j][i]) > fabs(a[lead][i]))
				lead = j;
		}

		if (lead != i) {
			double tv   = y[lead];
			y[lead]     = y[i];
			y[i]        = tv;

			double *tvm = a[i];
			a[i]        = a[lead];
			a[lead]     = tvm;
		}

#pragma omp parallel for schedule(dynamic, 16)
		for (int k = i + 1; k < N; ++k) {
			double fcr = a[k][i] / a[i][i];

			// #pragma omp simd
			for (int j = i; j < N; ++j) {
				a[k][j] -= fcr * a[i][j];
			}
			y[k] -= fcr * y[i];
		}
	}

	for (int i = N - 1; i >= 0; --i) {
		sol[i] = y[i];
		// #pragma omp simd
		for (int j = i + 1; j < N; ++j) {
			sol[i] -= a[i][j] * sol[j];
		}
		sol[i] /= a[i][i];
	}

	return sol;
}

long int
funcTask2(int **M, int N)
{
	double **L = calloc(N, sizeof(double *));
	double **U = calloc(N, sizeof(double *));

	for (int i = 0; i < N; ++i) {
		L[i] = calloc(N, sizeof(double));
		U[i] = calloc(N, sizeof(double));
	}

	for (int i = 0; i < N; ++i) {
		L[i][i] = 1;
	}

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < i; ++j) {
			double sum = 0;
			for (int k = 0; k < j; ++k) {
				sum += L[i][k] * U[k][j];
			}
			L[i][j] = (double)M[i][j] - sum;

			if (U[j][j] == 0) {
				printf("\ncould not do LU decompose, division by zero\n");
				return 0;
			}
			L[i][j] /= U[j][j];
		}

		for (int j = i; j < N; ++j) {
			double sum = 0;
			for (int k = 0; k < i; ++k) {
				sum += L[i][k] * U[k][j];
			}
			U[i][j] = (double)M[i][j] - sum;
		}
	}

	double res = 1;
	for (int i = 0; i < N; ++i) {
		res *= U[i][i];
	}

	for (int i = 0; i < N; ++i) {
		free(L[i]);
		free(U[i]);
	}
	free(L);
	free(U);

	return (long int)res;
}

#define COUNT_ARGS_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define COUNT_ARGS(...) \
COUNT_ARGS_IMPL(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

// matrix lib start
#define __matrix_destroy(...)                                             \
__matrix_destroy_impl(COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)
#define __matrix_ccreate(rows, ...)                                       \
	__matrix_ccreate_impl(rows, COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)

#define __matrix_val(m, i, j) ((m)->data[(i) * (m)->rows + (j)])

struct matrix {
	size_t  rows;
	double *data;
};

static void
__matrix_ccreate_impl(size_t rows, size_t count, ...)
{
	va_list args;
	va_start(args, count);

	for (size_t i = 0; i < count; ++i) {
		struct matrix *m = va_arg(args, struct matrix *);
		m->rows          = rows;

		m->data          = (double *)malloc(rows * rows * sizeof(double));
		memset(m->data, 0, rows * rows * sizeof(double));
	}
	va_end(args);
}

static void
__matrix_destroy_impl(size_t count, ...)
{
	va_list args;
	va_start(args, count);

	for (size_t i = 0; i < count; ++i) {
		struct matrix *m = va_arg(args, struct matrix *);
		m->rows          = 0;

		free(m->data);
		m->data = NULL;
	}
	va_end(args);
}
// matrix lib end

double **
funcTask3(float **A, float **B, int N)
{
	omp_set_num_threads(3);
	if (N <= 0)
		return NULL;

	size_t const  block      = 32;
	size_t        M          = ((N + block - 1) / block) * block;
	size_t        num_blocks = M / block;

	struct matrix a, b, c;
	__matrix_ccreate(M, &a, &b, &c);

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			__matrix_val(&a, i, j) = (double)A[i][j];
			__matrix_val(&b, i, j) = (double)B[i][j];
		}
	}

	double **res = malloc(sizeof(double *) * N);
	for (int i = 0; i < N; ++i)
		res[i] = malloc(sizeof(double) * N);

#pragma omp parallel for schedule(static)
	for (size_t bi = 0; bi < num_blocks; ++bi) {
		for (size_t bj = 0; bj < num_blocks; ++bj) {
			size_t row_base   = bi * block * M;
			size_t c_col_base = bj * block;

			for (size_t bk = 0; bk < num_blocks; ++bk) {
				size_t a_col_base = bk * block;
				size_t b_row_base = bk * block * M;
				size_t b_col_base = bj * block;

				for (size_t ii = 0; ii < block; ++ii) {
					double *a_ptr
					        = &a.data[row_base + ii * M + a_col_base];
					double *c_ptr
					        = &c.data[row_base + ii * M + c_col_base];

					__builtin_prefetch(a_ptr + 64, 0, 3);
					__builtin_prefetch(
					        &b.data[b_row_base + b_col_base], 0, 3);

#pragma unroll(8)
					for (size_t jj = 0; jj < block; ++jj) {
						double sum = 0.0;
#pragma unroll(8)
						for (size_t kk = 0; kk < block; ++kk) {
							sum += a_ptr[kk]
							       * b.data[b_row_base
							                + kk * M
							                + b_col_base
							                + jj];
						}
						c_ptr[jj] += sum;
					}
				}
			}
		}
	}

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			res[i][j] = __matrix_val(&c, i, j);
		}
	}
	__matrix_destroy(&a, &b, &c);

	return res;
}

const char *
funcLibInfoNickname()
{
	return "shgpavel";
}

const char *
funcLibInfoVersion()
{
	return "1.0";
}
