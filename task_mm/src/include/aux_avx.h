/* SPDX-License-Identifier: Apache-2.0 */

#ifndef AUX_AVX_H
#define AUX_AVX_H

#include <immintrin.h>

void avxreg_print(__m256d);

inline double avxreg_sum(__m256d r) {
	__m128d low = _mm256_castpd256_pd128(r);
	__m128d high = _mm256_extractf128_pd(r, 1);
	__m128d sum = _mm_add_pd(low, high);
	return _mm_cvtsd_f64(_mm_add_sd(sum, _mm_unpackhi_pd(sum, sum)));
}

inline double avxreg_sum512(__m512d r) {
	return _mm512_reduce_add_pd(r);
}

#endif
