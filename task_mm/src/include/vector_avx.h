/* SPDX-License-Identifier: Apache-2.0 */

#ifndef VECTOR_AVX_H
#define VECTOR_AVX_H

#include <immintrin.h>
#include "vector.h"

void   avxreg_print(__m256d); 
double vector_scalar_prod_avx(vector *, vector *);
double vector_norm2_avx(vector *, size_t);

inline double avxreg_sum(__m256d r) {
	__m256d sum = _mm256_hadd_pd(r, r);
	__m128d low  = _mm256_castpd256_pd128(sum);
	__m128d high = _mm256_extractf128_pd(sum, 1);
	__m128d sum2 = _mm_add_pd(low, high);
	return _mm_cvtsd_f64(sum2);
}

inline double avxreg_sum512(__m512d r) {
	return _mm512_reduce_add_pd(r);	
} 


#endif
