/* SPDX-License-Identifier: Apache-2.0 */

#include <jemalloc/jemalloc.h>
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

int compare_doubles(const void *a, const void *b) {
	double arg1 = *(const double *)a;
	double arg2 = *(const double *)b;

	if (arg1 < arg2) {
		return -1;
	} else if (arg1 > arg2) {
		return 1;
	}
	return 0;
}

void *aalloc(size_t size) {
	const size_t alignment = 64;
	if (size % alignment != 0) {
		size = (size + alignment - 1) & ~(alignment - 1);
	}
	return aligned_alloc(alignment, size);
}

const char outfp[] = "out.csv";
typedef void (*matrix_mult_func)(matrix *, matrix *, matrix *);

struct stats {
	double avg;
	double p95;
	double p1;
};

struct stats bench(matrix_mult_func multer, size_t n, size_t tests) {
	matrix a, b, c, btr;
	matrix_ccreate(&a, n, n, aalloc);
	matrix_ccreate(&b, n, n, aalloc);
	matrix_ccreate(&c, n, n, aalloc);
	matrix_ccreate(&btr, n, n, aalloc);

	generate_test(&a, &b, n);
	matrix_transpose(&btr, &b);

	vector times;
	vector_ccreate(&times, tests, aalloc);

	for (size_t i = 0; i < tests; ++i) {
		clock_t start = clock();
		if (multer == matrix_mult_fast) {
			multer(&a, &btr, &c);
		} else {
			multer(&a, &b, &c);
		}
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

	struct stats res = {avg, times.data[p95idx], low1p};
	matrix_destroy(&a, &b, &c, &btr);
	vector_destroy(&times);
	return res;
}

int main() {
	srand(time(NULL));

	/*
	FILE *csvout = fopen(outfp, "w");
	if (!csvout) {
		return -1;
	}

	int err = 0;
	err = fprintf(
	    csvout,
	    "n,Oavg,O95p,O1p,Navg,N95p,N1p,Ravg,R95p,R1p\n");
	if (err < 0) {
		return -1;
	}

	for (size_t n = 32; n < 1024; n += 32) {
		struct stats resO = bench(matrix_mult_openblas, n, 10);
		struct stats resN = bench(matrix_mult_naive, n, 10);
		struct stats resR = bench(matrix_mult_fast, n, 10);

		fprintf(csvout,
		        "%zu,%lg,%lg,%lg,%lg,%lg,%lg,%lg"
		        ",%lg,%lg\n",
		        n, resO.avg, resO.p95, resO.p1, resN.avg, resN.p95,
		        resN.p1, resR.avg, resR.p95, resR.p1);
	}

	err = fclose(csvout);
	if (err) {
		return -2;
	}

	return 0;
	*/

	size_t n;
	scanf("%zu", &n);
	
	matrix a, b, c, btr, c2;
	matrix_ccreate(&a, n, n, aalloc);
	matrix_ccreate(&b, n, n, aalloc);
	matrix_ccreate(&c, n, n, aalloc);
	matrix_ccreate(&c2, n, n, aalloc);

	matrix_ccreate(&btr, n, n, aalloc);

	for (size_t i = 0; i < a.rows; ++i) {
		for (size_t j = 0; j < a.cols; ++j) {
			scanf("%lf", &matrix_val(&a, i, j));
		}
	}

	for (size_t i = 0; i < b.rows; ++i) {
		for (size_t j = 0; j < b.cols; ++j) {
			scanf("%lf", &matrix_val(&b, i, j));
		}
	}

	matrix_print(&a);
	printf("\n");

	matrix_transpose(&btr, &b);
	matrix_print(&btr);
	printf("\n");

	matrix_mult_fast(&a, &btr, &c);
	matrix_print(&c);

	printf("\n");
	matrix_mult_openblas(&a, &b, &c2);
	matrix_print(&c2);

	matrix_destroy(&a, &b, &c, &btr, &c2);
}
