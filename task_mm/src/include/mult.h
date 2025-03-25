/* SPDX-License-Identifier: Apache-2.0 */

#ifndef MULT_H
#define MULT_H

#include "matrix.h"

void matrix_mult_mkl(matrix const *, matrix const *, matrix const *);

void matrix_mult_naive(matrix const *, matrix const *, matrix const *);

void matrix_mult_fast(matrix const *, matrix const *, matrix const *);

#endif
