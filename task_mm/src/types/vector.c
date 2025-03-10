/* SPDX-License-Identifier: Apache-2.0 */

#include <jemalloc/jemalloc.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <vector.h>

void vector_create(vector *v, size_t capacity, allocator a) {
	v->data = (double *)a(capacity * sizeof(double));
	v->size = 0;
	v->capacity = capacity;
}

void vector_ccreate(vector *v, size_t capacity, allocator a) {
	vector_create(v, capacity, a);
	memset(v->data, 0, capacity * sizeof(double));
	v->size = capacity;
}

void vector_create_copy(vector *dest, vector *src, allocator a) {
	vector_ccreate(dest, src->capacity, a);
	for (size_t i = 0; i < src->size; ++i) {
		dest->data[i] = src->data[i];
	}
	dest->size = src->size;
}

void vector_resize(vector *v, size_t new_capacity, reallocator rea) {
	v->capacity = new_capacity;
	v->data = rea(v->data, new_capacity * sizeof(double));
}

void vector_destroy_impl(size_t count, ...) {
	va_list args;
	va_start(args, count);

	for (size_t i = 0; i < count; ++i) {
		vector *v = va_arg(args, vector *);
		v->size = 0;
		v->capacity = 0;
		free(v->data);
		v->data = NULL;
	}
	va_end(args);
}

void vector_push(vector *v, double atad, reallocator rea) {
	if (v->size >= v->capacity) {
		vector_resize(v, (v->capacity + 1) + v->capacity / 2, rea);
	}

	v->data[v->size] = atad;
	v->size++;
}

void vector_assign(vector *v, vector *c, allocator a, reallocator rea) {
	if (v->size > c->size) {
		vector_destroy_impl(1, v);
		vector_create_copy(v, c, a);
	} else {
		for (size_t i = 0; i < v->size; ++i) {
			v->data[i] = c->data[i];
		}

		if (v->size < c->size) {
			for (size_t i = v->size; i < c->size; ++i) {
				vector_push(v, c->data[i], rea);
			}
		}
	}
}

void vector_print(vector *v) {
	if (v != NULL) {
		for (size_t i = 0; i < v->size; ++i) {
			printf("%g ", v->data[i]);
		}
		printf("\n");
	}
}

void vector_swap_eff(vector *v, vector *c) {
	size_t tmp_st = v->size;
	v->size = c->size;
	c->size = tmp_st;

	tmp_st = v->capacity;
	v->capacity = c->capacity;
	c->capacity = tmp_st;

	double *tmp = v->data;
	v->data = c->data;
	c->data = tmp;
}

void vector_fill_smth(vector *v, double x) {
	for (size_t i = 0; i < v->capacity; ++i) {
		v->data[i] = x;
	}
}
