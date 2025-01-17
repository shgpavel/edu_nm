/* SPDX-License-Identifier: Apache-2.0 */

#ifndef VECTOR_AVX_H
#define VECTOR_AVX_H

#include <immintrin.h>
#include "vector.h"

void print_avxreg(__m256d);
double vector_scalar_prod_avx(vector *, vector *);
double vector_norm2_avx(vector *, size_t);

#endif
