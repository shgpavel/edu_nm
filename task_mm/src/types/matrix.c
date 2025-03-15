/* SPDX-License-Identifier: Apache-2.0 */

#include "matrix.h"

#include <jemalloc/jemalloc.h>
#include <stdarg.h>
#include <stdio.h>

#include "vector.h"

void matrix_create(matrix *m, size_t rows, size_t cols, allocator a) {
	m->rows = rows;
	m->cols = cols;
	m->data = (vector *)a(sizeof(vector));
	vector_create(m->data, m->rows * m->cols, a);
}

void matrix_ccreate_impl(size_t rows, size_t cols, allocator alloc,
                         size_t count, ...) {
	va_list args;
	va_start(args, count);

	for (size_t i = 0; i < count; ++i) {
		matrix *m = va_arg(args, matrix *);
		m->rows = rows;
		m->cols = cols;

		m->data = (vector *)alloc(sizeof(vector));
		vector_ccreate(m->data, m->rows * m->cols, alloc);
	}
	va_end(args);
}

void matrix_create_copy(matrix *dest, matrix *src, allocator a) {
	if (src != NULL) {
		dest->rows = src->rows;
		dest->cols = src->cols;
		dest->data = (vector *)malloc(sizeof(vector));
		vector_create_copy(dest->data, src->data, a);
	}
}

void matrix_push(matrix *m, double adat, reallocator rea) {
	vector_push(m->data, adat, rea);
}

void matrix_destroy_impl(size_t count, ...) {
	va_list args;
	va_start(args, count);

	for (size_t i = 0; i < count; ++i) {
		matrix *m = va_arg(args, matrix *);
		m->rows = 0;
		m->cols = 0;

		vector_destroy(m->data);
		free(m->data);
		m->data = NULL;
	}
	va_end(args);
}

void matrix_print(matrix *m) {
	int max_len = 0;
	for (size_t i = 0; i < m->rows; ++i) {
		for (size_t j = 0; j < m->cols; ++j) {
			double val = matrix_val(m, i, j);
			int len = snprintf(NULL, 0, "%g", val);
			if (len > max_len) {
				max_len = len;
			}
		}
	}

	for (size_t i = 0; i < m->rows; ++i) {
		for (size_t j = 0; j < m->cols; ++j) {
			printf("%*lg ", max_len, matrix_val(m, i, j));
		}
		printf("\n");
	}
}

void matrix_transpose(matrix *dest, matrix *src) {
	if (src->cols != dest->rows || src->rows != dest->cols) {
		return;
	}

	for (size_t i = 0; i < dest->rows; ++i) {
		for (size_t j = 0; j < dest->cols; ++j) {
			matrix_val(dest, i, j) = matrix_val(src, j, i);
		}
	}
}
