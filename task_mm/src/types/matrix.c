/* SPDX-License-Identifier: Apache-2.0 */

#include <jemalloc/jemalloc.h>
#include <math.h>
#include <matrix.h>
#include <stdarg.h>
#include <stdio.h>
#include <vector.h>

void matrix_create(matrix *m, size_t rows, size_t cols, allocator a) {
	m->rows = rows;
	m->cols = cols;
	m->data = (vector *)a(sizeof(vector));
	vector_create(m->data, m->rows * m->cols, a);
}

void matrix_ccreate(matrix *m, size_t rows, size_t cols, allocator a) {
	m->rows = rows;
	m->cols = cols;
	m->data = (vector *)a(sizeof(vector));
	vector_ccreate(m->data, m->rows * m->cols, a);
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
	for (size_t i = 0; i < m->rows; ++i) {
		for (size_t j = 0; j < m->cols; ++j) {
			printf("%lg ", matrix_val(m, i, j));
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

int matrix_equal(matrix *a, matrix *b) {
	if (a->rows != b->rows || a->cols != b->cols) {
		return -1;
	}

	int rly = 1;
	for (size_t i = 0; i < a->rows; ++i) {
		for (size_t j = 0; j < a->cols; ++j) {
			if (fabs(matrix_val(a, i, j) - matrix_val(b, i, j)) >
			    1e-3) {
				rly = -1;
			}
		}
	}

	if (rly > 0) {
		return 1;
	}
	return -1;
}
