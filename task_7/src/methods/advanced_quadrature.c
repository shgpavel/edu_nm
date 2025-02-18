/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Contains advanced quadrature formulas with 3-point interpolation.
 * - Newton-Cotes
 * - Gauss
 */

#include <jemalloc/jemalloc.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "cube_eq.h"
#include "func.h"
#include "matrix.h"
#include "qr_avx.h"
#include "vector.h"
#include "vector_avx.h"

#define N 3

int combinations(int n, int k) {
	int Bm[n + 1][n + 1];
	for (int i = 0; i <= n; ++i) {
		Bm[i][0] = 1;
		Bm[i][i] = 1;
		for (int j = 1; j < i; ++j) {
			Bm[i][j] = Bm[i - 1][j - 1] + Bm[i - 1][j];
		}
	}
	return Bm[n][k];
}

double compute_binomial(double m, double p, double a, double t) {
	double total = 0.0;
	for (size_t i = 0; i < (size_t)m + 1; ++i) {
		int c = combinations((int)m, i);
		total +=
		    c * pow(a, i) * pow(t - a, m - i + p + 1) / (m - i + p + 1);
	}
	return total;
}

double calc_moment(size_t i, double x) {
	return compute_binomial((double)i, -ALPHA, A, x);
}

double newton_cotes(double a, double b) {
	vector points;
	vector_cinit(&points, N, sizeof(double));
	vector_val(&points, 0) = a;
	vector_val(&points, 1) = (a + b) / 2;
	vector_val(&points, 2) = b;

	vector fvals;
	vector_cinit(&fvals, N, sizeof(double));
	for (size_t i = 0; i < fvals.size; ++i) {
		vector_val(&fvals, i) = func(vector_val(&points, i));
	}

	vector moments;
	vector_cinit(&moments, N, sizeof(double));
	for (size_t i = 0; i < N; ++i) {
		vector_val(&moments, i) = calc_moment(i, b) - calc_moment(i, a);
	}

	matrix vtt;
	matrix_cinit(&vtt, N, N, sizeof(double));
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			matrix_val(&vtt, i, j) = pow(vector_val(&points, j), i);
		}
	}
	vector_free(&points);

	vector *sol = qr_avx(&vtt, &moments);
	matrix_free(&vtt);
	vector_free(&moments);

	double res = vector_scalar_prod_avx(sol, &fvals);

	vector_free(&fvals);
	vector_free(sol);
	free(sol);

	return res;
}

double gauss_quad(double a, double b) {
	vector moments;
	vector moments_rs;
	vector_cinit(&moments, 2 * N, sizeof(double));
	vector_cinit(&moments_rs, N, sizeof(double));

	for (size_t i = 0; i < 2 * N; ++i) {
		vector_val(&moments, i) = calc_moment(i, b) - calc_moment(i, a);
	}
	for (size_t i = 0; i < N; ++i) {
		vector_val(&moments_rs, i) = -vector_val(&moments, i + N);
	}
	// vector_print(&moments);
	// vector_print(&moments_rs);

	matrix mm;
	matrix_cinit(&mm, N, N, sizeof(double));
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			matrix_val(&mm, i, j) = vector_val(&moments, i + j);
		}
	}
	// matrix_print(&mm);

	vector *sol = qr_avx(&mm, &moments_rs);
	matrix_free(&mm);
	vector_free(&moments_rs);

	vector copoly;
	vector_cinit(&copoly, N + 1, sizeof(double));
	for (size_t i = 0; i < copoly.size - 1; ++i) {
		vector_val(&copoly, i) = vector_val(sol, i);
	}
	vector_val(&copoly, copoly.size - 1) = 1;
	vector_free(sol);
	free(sol);
	// vector_print(&copoly);

	vector roots = cube_eq_roots(&copoly);
	vector_free(&copoly);
	// vector_print(&roots);

	vector fvals;
	vector_cinit(&fvals, N, sizeof(double));
	for (size_t i = 0; i < fvals.size; ++i) {
		vector_val(&fvals, i) = func(vector_val(&roots, i));
	}

	matrix vtt;
	vector moments_re;
	matrix_cinit(&vtt, N, N, sizeof(double));
	vector_cinit(&moments_re, N, sizeof(double));

	for (size_t i = 0; i < N; ++i) {
		vector_val(&moments_re, i) = vector_val(&moments, i);
	}
	vector_free(&moments);

	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			matrix_val(&vtt, i, j) = pow(vector_val(&roots, j), i);
		}
	}

	vector *solf = qr_avx(&vtt, &moments_re);
	matrix_free(&vtt);
	vector_free(&moments_re);

	double res = vector_scalar_prod_avx(solf, &fvals);

	vector_free(&fvals);
	vector_free(&roots);
	vector_free(solf);
	free(solf);

	return res;
}

double composite_newton_cotes(size_t n) {
	double step = (B - A) / n;
	double total = 0.0;
	for (size_t i = 0; i < n; ++i) {
		double ai = A + i * step;
		double bi = ai + step;
		total += newton_cotes(ai, bi);
	}
	return total;
}

double composite_gauss(size_t n) {
	double step = (B - A) / n;
	double total = 0.0;
	for (size_t i = 0; i < n; ++i) {
		double ai = A + i * step;
		double bi = ai + step;
		total += gauss_quad(ai, bi);
	}
	return total;
}

int main() {
	printf("Newton_Cotes: %lf\n", composite_newton_cotes(10));
	printf("Gauss_Quad: %lf\n", composite_gauss(10));
	return 0;
}
