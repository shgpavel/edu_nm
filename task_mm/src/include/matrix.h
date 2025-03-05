/* SPDX-License-Identifier: Apache-2.0 */

#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include "vector.h"

typedef struct matrix_s {
  size_t rows;
  size_t cols;
  vector *data;
} matrix;

#define matrix_destroy(...) \
	matrix_destroy_impl(COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)
#define matrix_val(m, i, j) ((m)->data->data[(i) * (m)->cols + (j)])


void   matrix_create(matrix *, size_t, size_t, allocator);
void   matrix_ccreate(matrix *, size_t, size_t, allocator);
void   matrix_create_copy(matrix *, matrix *, allocator);

void   matrix_push(matrix *, double, reallocator);
void   matrix_destroy_impl(size_t, ...);
void   matrix_print(matrix *);

void   matrix_on_matrix(matrix *, matrix *, matrix *);
void   matrix_transpose(matrix *, matrix *);
int    matrix_equal(matrix *, matrix *);
/*
void matrix_change(matrix *, size_t, size_t, void *);
void matrix_swap(matrix *, size_t, size_t, size_t, size_t);
void matrix_row_swap(matrix *, size_t, size_t);
void matrix_delete(matrix *, size_t, size_t);
void matrix_copy(matrix *, matrix *);
void matrix_fill_smth(matrix *, double);
double matrix_norm_inf(matrix *);
vector *matrix_on_vector(matrix *, vector *);
void matrix_normalize_vect(matrix *, vector *);
matrix *matrix_on_matrix(matrix *, matrix *);
matrix *matrix_inverse(matrix *);
*/

#endif
