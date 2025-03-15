/* SPDX-License-Identifier: Apache-2.0 */

#ifndef MM_OB
#define MM_OB

#include "matrix.h"

void matrix_mult_mkl(matrix *, matrix *, matrix *);
void matrix_mult_naive(matrix *, matrix *, matrix *);
void matrix_mult_fast(matrix *, matrix *, matrix *);

#endif
