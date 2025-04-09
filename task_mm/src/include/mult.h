/* SPDX-License-Identifier: Apache-2.0 */

#ifndef MULT_H
#define MULT_H

#include "matrix.h"

void matrix_mult_mkl(matrix *, matrix *, matrix *);

void matrix_mult_naive(matrix *, matrix *, matrix *);

void matrix_mult_fast(matrix *, matrix *, matrix *);

#endif
