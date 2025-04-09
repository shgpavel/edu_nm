/* SPDX-License-Identifier: Apache-2.0 */

#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

typedef void *(*allocator)(size_t);
typedef void *(*reallocator)(void *, size_t);

#define COUNT_ARGS_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define COUNT_ARGS(...) \
	COUNT_ARGS_IMPL(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

#define matrix_destroy(...)																	\
	matrix_destroy_impl(COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)
#define matrix_ccreate(rows, cols, allocator, reallocator, ...)				\
	matrix_ccreate_impl(rows, cols, allocator, reallocator, COUNT_ARGS(__VA_ARGS__), \
	                    __VA_ARGS__)
#define matrix_val(m, i, j) ((m)->data[(i) * (m)->cols + (j)])
#define matrix_direct(m, i) ((m)->data[i])

typedef struct matrix_s {
	size_t rows;
	size_t cols;
	allocator alloc;
	reallocator realloc;
	double *data;
} matrix;

void matrix_ccreate_impl(size_t, size_t, allocator, reallocator, size_t, ...);
void matrix_destroy_impl(size_t, ...);
void matrix_print(matrix *);

void matrix_transpose(matrix *, matrix *);

void matrix_resize_specc(matrix *, size_t);
void matrix_resize_specr(matrix *, size_t);

#endif
