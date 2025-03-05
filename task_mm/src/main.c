/* SPDX-License-Identifier: Apache-2.0 */

#include <jemalloc/jemalloc.h>
#include <stdio.h>
#include <time.h>

#include "include/matrix.h"
#include "matrix.h"
#include "mm_ob.h"
#include "vector.h"
#include "vector_avx.h"

int generate_test(matrix *a, matrix *b, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			matrix_val(a, i, j) = rand();
			matrix_val(b, i, j) = (rand() - 97) / 1.2;
		}
	}
	return 0;
}

int compare_doubles(const void* a, const void* b) {
	double arg1 = *(const double*)a;
	double arg2 = *(const double*)b;
  
	if (arg1 < arg2) {
		return -1;
	} else if (arg1 > arg2) {
		return 1;
	}
	return 0;
}

typedef void (*matrix_mult_func)(matrix *, matrix *, matrix *);
int bench(matrix_mult_func multyk, size_t n, size_t tests) {
	matrix a, b, c;
	matrix_ccreate(&a, n, n, malloc);
	matrix_ccreate(&b, n, n, malloc);
	matrix_ccreate(&c, n, n, malloc);
	
	generate_test(&a, &b, n);

	vector times;
	vector_ccreate(&times, tests, malloc);
	
	for (size_t i = 0; i < tests; ++i) {
		clock_t start = clock();
		multyk(&a, &b, &c);
		clock_t end = clock();
		double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
		times.data[i] = time_spent;
	}

	double avg = times.data[0];
	for (size_t i = 1; i < times.size; ++i) {
		avg += times.data[i] / i - avg / i;
	}

	size_t low1p_cnt = (1 > tests / 100) ? 1 : tests / 100;
	qsort(times.data, times.size, sizeof(double), compare_doubles);
	double low1p = times.data[times.size - low1p_cnt - 1];
	for (size_t i = times.size - low1p_cnt; i < times.size; ++i) {
		low1p += times.data[i] / i - low1p / i;
	}

	size_t p95idx = (size_t)(0.95 * times.size + 0.5);
	if (p95idx == times.size) {
		--p95idx;
	}
	printf("times <sec>\tavg: %lg\t95p: %lg\t1p: %lg\n",
				 avg, times.data[p95idx], low1p);

	matrix c_et;
	matrix_ccreate(&c_et, n, n, malloc);
	matrix_mult_openblas(&a, &b, &c_et);
	
	if (matrix_equal(&c, &c_et) < 0) {
		printf("Error: mult function results are incorrect!\n");
	}
	
	matrix_destroy(&a, &b, &c, &c_et);
	vector_destroy(&times);
	return 0;
}


int main() {
	srand(time(NULL));

	printf(";; OpenBLAS res\n");
	bench(matrix_mult_openblas, 1000, 100);
	printf(";;\n\n");

	printf(";; MM naive res\n");
	bench(matrix_mult_naive, 50, 10);
	printf(";;\n\n");
	
	printf(";; MM transp res\n");
	bench(matrix_mult_fast, 50, 10);
	printf(";;\n");
	
	return 0;
}
