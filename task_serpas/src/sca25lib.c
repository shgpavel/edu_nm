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

/*
// gauss low(x86) high(e2k)
double *
funcTask1(double **a, double *y, int N)
{
        omp_set_num_threads(3);

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

#pragma omp parallel for schedule(static)
                for (int k = i + 1; k < N; ++k) {
                        double fcr = a[k][i] / a[i][i];

                        for (int j = i + 1; j < N; ++j) {
                                a[k][j] -= fcr * a[i][j];
                        }
                        y[k] -= fcr * y[i];
                }
        }

        for (int i = N - 1; i >= 0; --i) {
                double sum = 0.0;
                for (int j = i + 1; j < N; ++j) {
                        sum += a[i][j] * y[j];
                }
                y[i] = (y[i] - sum) / a[i][i];
        }

        return y;
}
*/

double *
funcTask1(double **a, double *y, int N)
{
	omp_set_num_threads(3);

	double **local_a = a;
	double  *local_y = y;

	for (int i = 0; i < N; ++i) {
		int    lead    = i;
		double max_val = fabs(local_a[i][i]);

		for (int j = i + 1; j < N; ++j) {
			double abs_val = fabs(local_a[j][i]);
			if (abs_val > max_val) {
				max_val = abs_val;
				lead    = j;
			}
		}

		if (lead != i) {
			double *tmp_row = local_a[i];
			local_a[i]      = local_a[lead];
			local_a[lead]   = tmp_row;

			double tmp_val  = local_y[i];
			local_y[i]      = local_y[lead];
			local_y[lead]   = tmp_val;
		}

		double inv_aii = 1.0 / local_a[i][i];

#pragma omp parallel for schedule(static)
		for (int k = i + 1; k < N; ++k) {
			double  factor  = local_a[k][i] * inv_aii;

			double *a_i_row = local_a[i] + i + 1;
			double *a_k_row = local_a[k] + i + 1;

			int     j       = i + 1;
			for (; j < N && ((uintptr_t)(a_k_row) % 64 != 0); ++j) {
				a_k_row[j - i - 1] -= factor * a_i_row[j - i - 1];
			}

			for (; j <= N - 4; j += 4) {
				a_k_row[j - i - 1] -= factor * a_i_row[j - i - 1];
				a_k_row[j - i]     -= factor * a_i_row[j - i];
				a_k_row[j - i + 1] -= factor * a_i_row[j - i + 1];
				a_k_row[j - i + 2] -= factor * a_i_row[j - i + 2];
			}

			for (; j < N; ++j) {
				a_k_row[j - i - 1] -= factor * a_i_row[j - i - 1];
			}

			local_y[k] -= factor * local_y[i];
		}
	}

	for (int i = N - 1; i >= 0; --i) {
		double  sum     = 0.0;
		double *a_i_row = local_a[i] + i + 1;
		double *y_ptr   = local_y + i + 1;

		for (int j = i + 1; j <= N - 4; j += 4) {
			sum += a_i_row[j - i - 1] * y_ptr[j - i - 1];
			sum += a_i_row[j - i] * y_ptr[j - i];
			sum += a_i_row[j - i + 1] * y_ptr[j - i + 1];
			sum += a_i_row[j - i + 2] * y_ptr[j - i + 2];
		}

		for (int j = (N - (N - i - 1) % 4); j < N; ++j) {
			sum += local_a[i][j] * local_y[j];
		}

		local_y[i] = (local_y[i] - sum) / local_a[i][i];
	}

	return local_y;
}

long int
funcTask2(int **M, int N)
{
	double *A = (double *)malloc(N * N * sizeof(double));

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
				printf("\ncould not do LU decompose, division by zero\n");
				free(A);
				return 0;
			}
			A[j * N + i] = (A[j * N + i] - sum) / A[i * N + i];
		}
	}

	double det = 1.0;
	for (int i = 0; i < N; ++i) {
		det *= A[i * N + i];
	}
	free(A);

	return (long int)det;
}

double **
funcTask3(float **A, float **B, int N)
{
	omp_set_num_threads(3);
	if (N <= 0)
		return NULL;

	const int block = 32;

	double  **res   = malloc(sizeof(double *) * N);
	if (!res)
		return NULL;

	double *data_block = calloc(N * N, sizeof(double));
	if (!data_block) {
		free(res);
		return NULL;
	}

	for (int i = 0; i < N; ++i) {
		res[i] = data_block + i * N;
	}

#pragma omp parallel for schedule(static) collapse(2)
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
