/* SPDX-License-Identifier: Apache-2.0 */

#ifndef MM_OB
#define MM_OB

#include "matrix.h"

void matrix_mult_mkl(matrix const* restrict, matrix const* restrict,
                     matrix const* restrict);

void matrix_mult_naive(matrix const* restrict, matrix const* restrict,
                       matrix const* restrict);

void matrix_mult_fast(matrix const* restrict, matrix const* restrict,
                      matrix const* restrict);

#endif
