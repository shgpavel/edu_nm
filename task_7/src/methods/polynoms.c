/* SPDX-License-Identifier: Apache-2.0 */
#include <jemalloc/jemalloc.h>

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>

#include "vector.h"
#include "polynoms.h"

void fft(complex double *x, size_t n, size_t inverse) {
	if (n <= 1) return;

	complex double even[n / 2];
	complex double odd[n / 2];
	for (size_t i = 0; i < n / 2; i++) {
		even[i] = x[2 * i];
		odd[i] = x[2 * i + 1];
	}

	fft(even, n / 2, inverse);
	fft(odd, n / 2, inverse);

	complex double W_N = cexp((2 * M_PI * I) / (double)n * (inverse ? 1 : -1));
	complex double W = 1;
	for (size_t k = 0; k < n / 2; k++) {
		complex double t = W * odd[k];
		x[k] = even[k] + t;
		x[k + n / 2] = even[k] - t;
		W *= W_N;
	}
}

vector poly_mult(size_t count, ...) {

	vector res;
	va_list ap;
	size_t total_size = 0;
	
	va_start(ap, count);
	for (size_t i = 0; i < count; ++i) {
		vector *poly = va_arg(ap, vector *);
		total_size += poly->size;
	}
	va_end(ap);
				
	size_t result_deg = total_size - (count - 1);
	size_t n = 1;
	while (n < result_deg) {
		n = n << 1;
	}

	vector_init(&res, result_deg, sizeof(double));

	
	complex double *fft_poly =
		(complex double *)calloc(n, sizeof(complex double));
	for (size_t i = 0; i < n; ++i) {
		fft_poly[i] = 1.0;
	}

	va_start(ap, count);
	for (size_t i = 0; i < count; ++i) {
		vector *poly = va_arg(ap, vector *);
		complex double *temp_poly =
			(complex double *)calloc(n, sizeof(complex double));
		for (size_t j = 0; j < poly->size; ++j) {
			temp_poly[j] = vector_val(poly, j);
		}
		fft(temp_poly, n, 0);
		for (size_t j = 0; j < n; j++) {
			fft_poly[j] *= temp_poly[j];
		}
		free(temp_poly);
	}
	va_end(ap);
	
  fft(fft_poly, n, 1);

	for (size_t i = 0; i < result_deg; ++i) {
		double tmp = (creal(fft_poly[i]) / (double)n);
		vector_push(&res, &tmp);
	}
	free(fft_poly);
	return res;
}

vector poly_sum(size_t count, ...) {

	va_list args, args_copy;
	va_copy(args_copy, args);
	
	size_t max_size = 0;
	va_start(args_copy, count);
	for (size_t i = 0; i < count; ++i) {
		vector *vec = va_arg(args_copy, vector *);
		max_size = vec->size > max_size ? vec->size : max_size;
	}
	va_end(args_copy);

	vector res;
	vector_init(&res, max_size, sizeof(double));
	vector_fill_smth(&res, 0.0);
	
	va_start(args, count);
	for (size_t i = 0; i < count; ++i) {
		vector *vec = va_arg(args, vector *);
		for (size_t j = 0; j < vec->size; ++j) {
			vector_val(&res, j) += vector_val(vec, j);
		}
	}
	va_end(args);

	return res;
}

double poly_val(vector *v, double point) {
  double res = 0.0;
  for (size_t i = 0; i < v->size; ++i) {
    res += vector_val(v, i) * pow(point, (double)i);
  }
  return res;
}
