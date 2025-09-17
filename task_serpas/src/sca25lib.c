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

#include <immintrin.h>

#define COUNT_ARGS_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define COUNT_ARGS(...) \
COUNT_ARGS_IMPL(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

// vector lib start
#define __vector_ccreate(capacity, ...)                                       \
	__vector_ccreate_impl(capacity, COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)
#define __vector_destroy(...) __vector_destroy_impl(COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)

#define __vector_val(v, i)    ((v)->data[(i)])

struct vector {
	size_t  size;
	size_t  capacity;
	double *data;
};

static void
__vector_ccreate_impl(size_t capacity, size_t count, ...)
{
	va_list args;
	va_start(args, count);

	for (size_t i = 0; i < count; ++i) {
		struct vector *v = va_arg(args, struct vector *);
		v->size          = capacity;
		v->data          = (double *)malloc(capacity * sizeof(double));
		memset(v->data, 0, capacity * sizeof(double));
		v->capacity = capacity;
	}
	va_end(args);
}

static void
__vector_destroy_impl(size_t count, ...)
{
	va_list args;
	va_start(args, count);

	for (size_t i = 0; i < count; ++i) {
		struct vector *v = va_arg(args, struct vector *);
		v->size          = 0;
		free(v->data);
		v->capacity = 0;
		v->data     = NULL;
	}
	va_end(args);
}

static void
__vector_swap_el(struct vector *v, size_t i, size_t k)
{
	if (i < v->size && k < v->size) {
		double tv  = v->data[k];
		v->data[k] = v->data[i];
		v->data[i] = tv;
	}
}

static double
__vector_norm2(struct vector *v, size_t k)
{
	double res = 0.0;
	// #pragma omp parallel for reduction(+ : res)
	for (size_t i = k; i < v->size; ++i) {
		res += v->data[i] * v->data[i];
	}
	return sqrt(res);
}
// vector lib end

// matrix lib start
#define __matrix_destroy(...) __matrix_destroy_impl(COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)
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

static void
__matrix_row_swap(struct matrix *m, size_t i, size_t j)
{
	if (j == i)
		return;

	double *tv = (double *)malloc(m->rows * sizeof(double));

	double *fo = m->data + i * m->rows;
	double *go = m->data + j * m->rows;

	memcpy(tv, go, m->rows * sizeof(double));
	memmove(go, fo, m->rows * sizeof(double));
	memcpy(fo, tv, m->rows * sizeof(double));

	free(tv);
}
// matrix lib end

// gauss low(x86) high(e2k)
double *
funcTask1(double **a, double *y, int N)
{
	omp_set_num_threads(3);

	struct matrix a_copy;
	struct vector y_copy, sol;
	__matrix_ccreate(N, &a_copy);
	__vector_ccreate(N, &y_copy, &sol);

	for (size_t i = 0; i < a_copy.rows; ++i) {
		for (size_t j = 0; j < a_copy.rows; ++j) {
			__matrix_val(&a_copy, i, j) = a[i][j];
		}
		__vector_val(&y_copy, i) = y[i];
	}

	for (size_t i = 0; i < a_copy.rows; ++i) {
		size_t lead = i;

		for (size_t j = i + 1; j < a_copy.rows; ++j) {
			if (fabs(__matrix_val(&a_copy, j, i))
			    > fabs(__matrix_val(&a_copy, lead, i)))
				lead = j;
		}

		if (lead != i) {
			__vector_swap_el(&y_copy, i, lead);
			__matrix_row_swap(&a_copy, i, lead);
		}

#pragma omp parallel for schedule(static)
		for (size_t k = i + 1; k < a_copy.rows; ++k) {
			double fcr = __matrix_val(&a_copy, k, i)
			             / __matrix_val(&a_copy, i, i);

			for (size_t j = i; j < a_copy.rows; ++j) {
				__matrix_val(&a_copy, k, j)
				        -= fcr * __matrix_val(&a_copy, i, j);
			}
			__vector_val(&y_copy, k) -= fcr * __vector_val(&y_copy, i);
		}
	}

	for (size_t i = a_copy.rows; i > 0;) {
		--i;
		__vector_val(&sol, i) = __vector_val(&y_copy, i);

		for (size_t j = i + 1; j < a_copy.rows; ++j) {
			__vector_val(&sol, i)
			        -= __matrix_val(&a_copy, i, j) * __vector_val(&sol, j);
		}
		__vector_val(&sol, i) /= __matrix_val(&a_copy, i, i);
	}

	__matrix_destroy(&a_copy);
	__vector_destroy(&y_copy);

	return sol.data;
}

/*
// gauss best(x86)
double *
funcTask1(double **a, double *y, int N)
{
        omp_set_num_threads(3);

        struct matrix a_copy;
        struct vector y_copy, sol;
        __matrix_ccreate(N, &a_copy);
        __vector_ccreate(N, &y_copy, &sol);

#pragma omp parallel for
        for (size_t i = 0; i < a_copy.rows; ++i) {
                for (size_t j = 0; j < a_copy.rows; ++j) {
                        __matrix_val(&a_copy, i, j) = a[i][j];
                }
                __vector_val(&y_copy, i) = y[i];
        }

        for (size_t i = 0; i < a_copy.rows; ++i) {
                size_t lead    = i;
                double max_val = fabs(__matrix_val(&a_copy, i, i));

#pragma omp parallel for reduction(max : max_val)
                for (size_t j = i + 1; j < a_copy.rows; ++j) {
                        double val = fabs(__matrix_val(&a_copy, j, i));
                        if (val > max_val) {
                                max_val = val;
                                lead    = j;
                        }
                }

                if (lead != i) {
                        __vector_swap_el(&y_copy, i, lead);
                        __matrix_row_swap(&a_copy, i, lead);
                }

                if (max_val < 1e-10) {
                        continue;
                }

                double diag = __matrix_val(&a_copy, i, i);

#pragma omp parallel for
                for (size_t k = i + 1; k < a_copy.rows; ++k) {
                        double fcr = __matrix_val(&a_copy, k, i) / diag;

#pragma omp simd
                        for (size_t j = i; j < a_copy.rows; ++j) {
                                __matrix_val(&a_copy, k, j)
                                        -= fcr * __matrix_val(&a_copy, i, j);
                        }
                        __vector_val(&y_copy, k) -= fcr * __vector_val(&y_copy,
i);
                }
        }

        for (size_t i = a_copy.rows; i > 0;) {
                --i;
                __vector_val(&sol, i) = __vector_val(&y_copy, i);

#pragma omp simd
                for (size_t j = i + 1; j < a_copy.rows; ++j) {
                        __vector_val(&sol, i) -= __matrix_val(&a_copy, i, j)
                                                 * __vector_val(&sol, j);
                }
                __vector_val(&sol, i) /= __matrix_val(&a_copy, i, i);
        }

        __matrix_destroy(&a_copy);
        __vector_destroy(&y_copy);

        return sol.data;
}
*/

// QR
double *
funcTask1_QR(double **a, double *y, int N)
{
	omp_set_num_threads(3);

	struct matrix a_copy;
	struct vector y_copy, sol, na, cacheV;
	__matrix_ccreate(N, &a_copy);
	__vector_ccreate(N, &y_copy, &na, &sol, &cacheV);

#pragma omp      parallel for schedule(static)
        for (size_t i = 0; i < a_copy.rows; ++i) {
                for (size_t j = 0; j < a_copy.rows; ++j) {
                        __matrix_val(&a_copy, i, j) = a[i][j];
                }
                __vector_val(&y_copy, i) = y[i];
        }

        for (size_t i = 0; i < a_copy.rows - 1; ++i) {
                for (size_t j = i; j < a_copy.rows; ++j) {
                        __vector_val(&na, j) = __matrix_val(&a_copy, j, i);
                }

                double norm = __vector_norm2(&na, i);
                if (norm == 0.0) {
                        continue;
                }

                double sign           = (__vector_val(&na, i) >= 0.0 ? 1.0 : -1.0);
                __vector_val(&na, i) += (sign * norm);

                double na_norm_2      = __vector_norm2(&na, i);
                if (na_norm_2 == 0.0) {
                        continue;
                }

                for (size_t j = i; j < a_copy.rows; ++j) {
                        __vector_val(&na, j) /= na_norm_2;
                }

#pragma omp parallel for schedule(static)
                for (size_t k = i; k < a_copy.rows; ++k) {
                        double sum = 0.0;
#pragma omp simd reduction(+ : sum) schedule(static)
                        for (size_t l = i; l < a_copy.rows; ++l) {
                                sum += __vector_val(&na, l) * __matrix_val(&a_copy, l, k);
                        }
                        __vector_val(&cacheV, k) = sum;
		     }

#pragma omp      parallel for schedule(static)
                for (size_t j = i; j < a_copy.rows; ++j) {
                        for (size_t k = i; k < a_copy.rows; ++k) {
                                __matrix_val(&a_copy, j, k)
                                        -= (2.0
                                            * __vector_val(&na, j)
                                            * __vector_val(&cacheV, k));
                        }
                }

                double sum_y = 0.0;
     #pragma omp parallel for reduction(+ : sum_y) schedule(static)
                for (size_t j = i; j < a_copy.rows; ++j) {
                        sum_y += __vector_val(&na, j) * __vector_val(&y_copy, j);
                }

                for (size_t j = i; j < a_copy.rows; ++j) {
                        __vector_val(&y_copy, j) -= (2.0 * __vector_val(&na, j) * sum_y);
                }
	     }

	     for (size_t i = a_copy.rows; i > 0;) {
                --i;
                __vector_val(&sol, i) = __vector_val(&y_copy, i);
                for (size_t j = a_copy.rows; j > i + 1;) {
                        --j;
                        __vector_val(&sol, i)
                                = __vector_val(&sol, i)
                                  - (__matrix_val(&a_copy, i, j) * __vector_val(&sol, j));
                }
                __vector_val(&sol, i)
                        = (__vector_val(&sol, i) / __matrix_val(&a_copy, i, i));
	     }

	     __matrix_destroy(&a_copy);
	     __vector_destroy(&y_copy, &na, &cacheV);

	     return sol.data;
}

long int
funcTask2(int **M, int N)
{
	// omp_set_num_threads(3);
	long **L = calloc(N, sizeof(long *));
	long **U = calloc(N, sizeof(long *));

	for (size_t i = 0; i < (size_t)N; ++i) {
		L[i] = calloc(N, sizeof(long));
		U[i] = calloc(N, sizeof(long));
	}

	for (size_t i = 0; i < (size_t)N; ++i) {
		L[i][i] = 1;
	}

	for (size_t i = 0; i < (size_t)N; ++i) {
		for (size_t j = 0; j < i; ++j) {
			long sum = 0;
			// #pragma omp parallel for reduction(+ : sum)
			for (size_t k = 0; k < j; ++k) {
				sum += L[i][k] * U[k][j];
			}
			L[i][j] = M[i][j] - sum;
			if (U[j][j] == 0) {
				printf("\ncould not do LU decompose, division by zero\n");
				return 0;
			}
			L[i][j] /= U[j][j];
		}

		// #pragma omp parallel for
		for (size_t j = i; j < (size_t)N; ++j) {
			long sum = 0;

			// #pragma omp parallel for reduction(+ : sum)
			for (size_t k = 0; k < i; ++k) {
				sum += L[i][k] * U[k][j];
			}
			U[i][j] = M[i][j] - sum;
		}
	}

	long res = 1;

	// #pragma omp parallel for reduction(* : res_u)
	for (size_t i = 0; i < (size_t)N; ++i) {
		res *= U[i][i];
	}

	for (size_t i = 0; i < (size_t)N; ++i) {
		free(L[i]);
		free(U[i]);
	}
	free(L);
	free(U);

	return res;
}

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

	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			__matrix_val(&a, i, j) = (double)A[i][j];
			__matrix_val(&b, i, j) = (double)B[i][j];
		}
	}

	double **res = malloc(sizeof(double *) * N);
	for (size_t i = 0; i < N; ++i)
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

	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
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
