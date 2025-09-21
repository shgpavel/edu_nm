/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Copyright (C) 2025 Pavel Shago <pavel@shago.dev>
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define _ISOC23_SOURCE

double *
funcTask1(double *restrict *restrict a, double *restrict y, int N)
{
	double tol   = 1e-3;
	int    maxit = 1000;

	omp_set_num_threads(3);

	double *x  = aligned_alloc(64, N * sizeof(double));
	double *r  = aligned_alloc(64, N * sizeof(double));
	double *r0 = aligned_alloc(64, N * sizeof(double));
	double *p  = aligned_alloc(64, N * sizeof(double));
	double *v  = aligned_alloc(64, N * sizeof(double));
	double *s  = aligned_alloc(64, N * sizeof(double));
	double *t  = aligned_alloc(64, N * sizeof(double));

	memset(x, 0, N * sizeof(double));
	memset(p, 0, N * sizeof(double));
	memset(v, 0, N * sizeof(double));

#pragma omp parallel for schedule(static)
	for (int i = 0; i < N; ++i) {
		double sum = 0.0;
		for (int j = 0; j < N; ++j) {
			sum += a[i][j] * x[j];
		}
		r[i]  = y[i] - sum;
		r0[i] = r[i];
	}

	double rho = 1.0, alpha = 1.0, omega = 1.0;

	for (int it = 0; it < maxit; ++it) {
		double rho_new = 0.0;

#pragma omp parallel for reduction(+ : rho_new) schedule(static)
		for (int i = 0; i < N; ++i) {
			rho_new += r0[i] * r[i];
		}

		if (fabs(rho_new) < 1e-20)
			break;

		double beta = (rho_new / rho) * (alpha / omega);

#pragma omp parallel for schedule(static)
		for (int i = 0; i < N; ++i) {
			p[i] = r[i] + beta * (p[i] - omega * v[i]);
		}

#pragma omp parallel for schedule(static)
		for (int i = 0; i < N; ++i) {
			double sum = 0.0;
			for (int j = 0; j < N; ++j) {
				sum += a[i][j] * p[j];
			}
			v[i] = sum;
		}

		double denom = 0.0;

#pragma omp parallel for reduction(+ : denom) schedule(static)
		for (int i = 0; i < N; ++i) {
			denom += r0[i] * v[i];
		}
		alpha = rho_new / denom;

#pragma omp parallel for schedule(static)
		for (int i = 0; i < N; ++i) {
			s[i] = r[i] - alpha * v[i];
		}

#pragma omp parallel for schedule(static)
		for (int i = 0; i < N; ++i) {
			double sum = 0.0;
			for (int j = 0; j < N; ++j) {
				sum += a[i][j] * s[j];
			}
			t[i] = sum;
		}

		double tt = 0.0, ts = 0.0;
#pragma omp parallel for reduction(+ : tt, ts) schedule(static)
		for (int i = 0; i < N; ++i) {
			tt += t[i] * t[i];
			ts += t[i] * s[i];
		}

		if (fabs(tt) < 1e-20)
			tt = 1e-20;
		omega = ts / tt;

#pragma omp parallel for schedule(static)
		for (int i = 0; i < N; ++i) {
			x[i] += alpha * p[i] + omega * s[i];
			r[i]  = s[i] - omega * t[i];
		}

		double norm = 0.0;
#pragma omp parallel for reduction(+ : norm) schedule(static)
		for (int i = 0; i < N; ++i) {
			norm += r[i] * r[i];
		}
		norm = sqrt(norm);

		if (norm < tol)
			break;

		rho = rho_new;
	}

	free(r);
	free(r0);
	free(p);
	free(v);
	free(s);
	free(t);

	return x;
}

long int
funcTask2(int *restrict *restrict M, int N)
{
	double A[N * N];

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			A[i * N + j] = (double)M[i][j];
		}
	}

	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j) {
			double sum = 0;
			for (int k = 0; k < i; ++k) {
				sum += A[i * N + k] * A[k * N + j];
			}
			A[i * N + j] -= sum;
		}

		for (int j = i + 1; j < N; ++j) {
			double sum = 0;
			for (int k = 0; k < i; ++k) {
				sum += A[j * N + k] * A[k * N + i];
			}
			if (A[i * N + i] == 0) {
				return 0;
			}
			A[j * N + i] = (A[j * N + i] - sum) / A[i * N + i];
		}
	}

	double det = 1.0;
	for (int i = 0; i < N; ++i) {
		det *= A[i * N + i];
	}

	return (long int)det;
}

double **
funcTask3(float *restrict *restrict A, float *restrict *restrict B, int N)
{
	omp_set_num_threads(3);
	if (N <= 0)
		return NULL;

	int const block                = 256;

	double *restrict *restrict res = aligned_alloc(64, sizeof(double *) * N);
	if (!res)
		return NULL;

	double *restrict data_block = aligned_alloc(64, N * N * sizeof(double));
	memset(data_block, 0, N * N * sizeof(double));
	if (!data_block) {
		free(res);
		return NULL;
	}

	for (int i = 0; i < N; ++i) {
		res[i] = data_block + i * N;
	}

#pragma omp parallel for schedule(static)
	for (int bi = 0; bi < N; bi += block) {
		for (int bj = 0; bj < N; bj += block) {
			int i_end = (bi + block) < N ? (bi + block) : N;
			int j_end = (bj + block) < N ? (bj + block) : N;

			for (int bk = 0; bk < N; bk += block) {
				int k_end = (bk + block) < N ? (bk + block) : N;

				for (int i = bi; i < i_end; ++i) {
					for (int k = bk; k < k_end; ++k) {
						double  a_val   = A[i][k];
						double *res_row = res[i] + bj;
						float  *b_row   = B[k] + bj;

						for (int j = 0; j < j_end - bj;
						     ++j) {
							res_row[j]
							        += a_val * b_row[j];
						}
					}
				}
			}
		}
	}

	return res;
}

#define INSERTION_THRESHOLD 32
#define PARALLEL_THRESHOLD  100000
#define CACHE_LINE_SIZE     32

static inline void
insertion_sort_desc(double *restrict arr, int left, int right)
{
	for (int i = left + 1; i <= right; i++) {
		double key = arr[i];
		int    j   = i - 1;

		while (j >= left + 3
		       && arr[j] < key
		       && arr[j - 1] < key
		       && arr[j - 2] < key
		       && arr[j - 3] < key) {
			arr[j + 1]  = arr[j];
			arr[j]      = arr[j - 1];
			arr[j - 1]  = arr[j - 2];
			arr[j - 2]  = arr[j - 3];
			j          -= 4;
		}

		while (j >= left && arr[j] < key) {
			arr[j + 1] = arr[j];
			j--;
		}
		arr[j + 1] = key;
	}
}

static void
merge_desc(double *restrict src, double *restrict dst, int left, int mid, int right)
{
	register double *p1   = &src[left];
	register double *p2   = &src[mid + 1];
	register double *end1 = &src[mid];
	register double *end2 = &src[right];
	register double *dest = &dst[left];

	__builtin_prefetch(p1, 0, 0);
	__builtin_prefetch(p2, 0, 0);

	while (p1 <= end1 && p2 <= end2) {
		register double v1 = *p1;
		register double v2 = *p2;

		if (v1 >= v2) {
			*dest++ = v1;
			p1++;
			__builtin_prefetch(p1 + 8, 0, 0);
		} else {
			*dest++ = v2;
			p2++;
			__builtin_prefetch(p2 + 8, 0, 0);
		}
	}

	while (p1 <= end1) {
		*dest++ = *p1++;
		__builtin_prefetch(p1 + 8, 0, 0);
	}

	while (p2 <= end2) {
		*dest++ = *p2++;
		__builtin_prefetch(p2 + 8, 0, 0);
	}
}

static void
parallel_merge_sort(double *restrict dst, double *restrict src, int left, int right,
                    int depth, double *temp)
{
	if (right - left <= INSERTION_THRESHOLD) {
		__builtin_prefetch(&src[left], 0, 0);
		__builtin_prefetch(&src[(left + right) / 2], 0, 0);
		__builtin_prefetch(&src[right], 0, 0);

		if (dst != src) {
			for (int i = left; i <= right; i++) {
				dst[i] = src[i];
			}
		}
		insertion_sort_desc(dst, left, right);
		return;
	}

	int     mid      = left + (right - left) / 2;
	double *next_dst = (dst == src) ? temp : src;

	__builtin_prefetch(&src[left], 0, 0);
	__builtin_prefetch(&src[mid], 0, 0);
	__builtin_prefetch(&src[right], 0, 0);

	if (depth > 0 && (right - left) > PARALLEL_THRESHOLD) {
#pragma omp task shared(src, dst, temp) firstprivate(left, mid, depth, next_dst)
		parallel_merge_sort(next_dst, src, left, mid, depth - 1, temp);

#pragma omp task shared(src, dst, temp) firstprivate(mid, right, depth, next_dst)
		parallel_merge_sort(next_dst, src, mid + 1, right, depth - 1, temp);

#pragma omp taskwait
	} else {
		parallel_merge_sort(next_dst, src, left, mid, 0, temp);
		parallel_merge_sort(next_dst, src, mid + 1, right, 0, temp);
	}

	merge_desc(next_dst, dst, left, mid, right);
}

void
parallel_sort(double *restrict arr, int n)
{
	double *temp  = aligned_alloc(CACHE_LINE_SIZE, n * sizeof(double));

	int     depth = 2;

#pragma omp parallel num_threads(3)
	{
#pragma omp single
		{
			parallel_merge_sort(arr, arr, 0, n - 1, depth, temp);
		}
	}

	free(temp);
}

void
funcTask4(double *restrict SortMe, int N)
{
	parallel_sort(SortMe, N);
}

const char *
funcLibInfoNickname()
{
	return "shgpavel";
}

const char *
funcLibInfoVersion()
{
	return "2.0";
}
