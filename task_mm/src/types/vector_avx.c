/* SPDX-License-Identifier: Apache-2.0 */

#include "vector_avx.h"

#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vector.h"

void avxreg_print(__m256d r) {
	double unpcd[4] __attribute__((aligned(64)));
	_mm256_store_pd(unpcd, r);
	for (size_t i = 0; i < 3; ++i) {
		printf("%lg ", unpcd[i]);
	}
	printf("%lg\n", unpcd[3]);
}

double vector_scalar_prod_avx(vector *v, vector *c) {
	size_t i = 0;
	double res = 0.0;
	if (v == NULL || v->data == NULL || c == NULL || c->data == NULL ||
	    c->size == 0 || v->size == 0) {
		return 0.0;
	}

	size_t n = v->size < c->size ? v->size : c->size;
	if (v->size < 4 || c->size < 4) {
		goto base;
	}

	__m256d sum = _mm256_setzero_pd();

	for (; i <= n - 4; i += 4) {
		double l1[4] __attribute__((aligned(64)));
		double l2[4] __attribute__((aligned(64)));

		for (size_t j = 0; j < 4; ++j) {
			l1[j] = v->data[i + j];
			l2[j] = c->data[i + j];
		}

		__m256d vec1 = _mm256_load_pd(l1);
		__m256d vec2 = _mm256_load_pd(l2);

		sum = _mm256_add_pd(sum, _mm256_mul_pd(vec1, vec2));
	}

	__m128d low = _mm256_castpd256_pd128(sum);
	__m128d high = _mm256_extractf128_pd(sum, 1);
	__m128d sum_r = _mm_add_pd(low, high);

	res =
	    _mm_cvtsd_f64(sum_r) + _mm_cvtsd_f64(_mm_unpackhi_pd(sum_r, sum_r));

base:
	for (; i < n; ++i) {
		res += v->data[i] * c->data[i];
	}
	return res;
}

double vector_norm2_avx(vector *v, size_t i) {
	double res = 0.0;
	if (v == NULL || v->data == NULL || v->size == 0) {
		return 0.0;
	}

	if (v->size < 4) {
		goto base;
	}

	__m256d sum = _mm256_setzero_pd();
	for (; i <= v->size - 4; i += 4) {
		double l1[4] __attribute__((aligned(64)));
		for (size_t j = 0; j < 4; ++j) {
			l1[j] = v->data[i + j];
		}

		__m256d vec1 = _mm256_load_pd(l1);
		sum = _mm256_add_pd(sum, _mm256_mul_pd(vec1, vec1));
	}

	__m128d low = _mm256_castpd256_pd128(sum);
	__m128d high = _mm256_extractf128_pd(sum, 1);
	__m128d sum_r = _mm_add_pd(low, high);

	res =
	    _mm_cvtsd_f64(sum_r) + _mm_cvtsd_f64(_mm_unpackhi_pd(sum_r, sum_r));

base:
	for (; i < v->size; ++i) {
		res += v->data[i] * v->data[i];
	}
	return sqrt(res);
}
