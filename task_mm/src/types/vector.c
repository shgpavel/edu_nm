/* SPDX-License-Identifier: Apache-2.0 */

#include "vector.h"

#include <jemalloc/jemalloc.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

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

/*
void vector_normalize(vector *v) {
        double norm = sqrt(vector_sclr_prod(v, v));
        for (size_t i = 0; i < v->size; ++i) {
                vector_val(v, i) /= norm;
        }
}

void vector_delete(vector *v, size_t idx) {
        if (index < v->size - 1) {
                free(vector_get(v, index));
                for (size_t i = index; i < v->size; ++i) {
                        v->data[i] = v->data[i + 1];
                }
                v->size--;
        } else if (index == v->size - 1) {
                free(vector_get(v, index));
                v->size--;
        }
}
*/

/*
void vector_swap(vector *v, size_t i, size_t j) {
        if (i < v->size && j < v->size) {
                void *tmp = (void *)malloc(v->type_size);
                memcpy(tmp, v->data[i], v->type_size);
                memcpy(v->data[i], v->data[j], v->type_size);
                memcpy(v->data[j], tmp, v->type_size);
                free(tmp);
        }
}

void vector_reverse(vector *v) {
        size_t start = 0, end = v->size - 1;
        double tmp;

        while (start < end) {
                tmp = vector_val(v, start);
                vector_val(v, start) = vector_val(v, end);
                vector_val(v, end) = tmp;

                ++start;
                --end;
        }
}

void vector_mult(vector *v, double a) {
        if (v->type_size == sizeof(double)) {
                for (size_t i = 0; i < v->size; ++i) {
                        vector_val(v, i) *= a;
                }
        }
}

double vector_diff(vector *x, vector *y) {
        if ((x->size == y->size) && (y->type_size == sizeof(double)) &&
            (x->type_size == sizeof(double))) {
                double result = 0.0;
                for (size_t i = 0; i < x->size; ++i) {
                        result += ((vector_val(x, i) - vector_val(y, i)) *
                                   (vector_val(x, i) - vector_val(y, i)));
                }
                return sqrt(result);
        }
        return 0.0;
}

double vector_sclr_prod(vector *v, vector *c) {
        double result = 0.0;
        for (size_t i = 0; i < v->size; ++i) {
                result += vector_val(v, i) * vector_val(c, i);
        }
        return result;
}

}
*/
