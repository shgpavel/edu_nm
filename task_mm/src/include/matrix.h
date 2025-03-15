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
#define matrix_ccreate(rows, cols, allocator, ...)		\
	matrix_ccreate_impl(rows, cols, allocator, COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)
#define matrix_val(m, i, j) ((m)->data->data[(i) * (m)->cols + (j)])
#define matrix_direct(m, i) ((m)->data->data[i])

void   matrix_create(matrix *, size_t, size_t, allocator);
void   matrix_ccreate_impl(size_t, size_t, allocator, size_t, ...);
void   matrix_create_copy(matrix *, matrix *, allocator);

void   matrix_push(matrix *, double, reallocator);
void   matrix_destroy_impl(size_t, ...);
void   matrix_print(matrix *);

void   matrix_transpose(matrix *, matrix *);

#endif
