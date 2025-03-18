/* SPDX-License-Identifier: Apache-2.0 */

#ifndef VECTOR_H
#define VECTOR_H

#include <stddef.h>

typedef struct vector_s {
	size_t size;
	size_t capacity;
	double *data;
} vector;

typedef void *(*allocator)(size_t);
typedef void *(*reallocator)(void *, size_t);

#define COUNT_ARGS_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define COUNT_ARGS(...) \
	COUNT_ARGS_IMPL(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#define vector_destroy(...) \
	vector_destroy_impl(COUNT_ARGS(__VA_ARGS__), __VA_ARGS__)

void vector_create(vector *, size_t, allocator);
void vector_ccreate(vector *, size_t, allocator);
void vector_create_copy(vector *, vector *, allocator);
void vector_print(vector *);

void vector_push(vector *, double, reallocator);
void vector_assign(vector *, vector *, allocator, reallocator);
void vector_destroy_impl(size_t count, ...);

#endif
