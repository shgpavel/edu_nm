/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Matrix multiplication
 */

#include "matrix.h"

void matrix_mult_fast(matrix *a, matrix *b, matrix *c) {
	/* TODO more checks */
	if (a->cols != b->rows) {
		return ;
	}

	for (size_t i = 0; i < a->rows; ++i) {
		for (size_t j = 0; j < b->cols; ++j) {
			double sum = 0.0;
			for (size_t k = 0; k < a->cols; ++k) {
				sum +=
				    matrix_val(a, i, k) * matrix_val(b, j, k);
			}
			matrix_val(c, i, j) = sum;
		}
	}	
}
