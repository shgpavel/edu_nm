/* SPDX-License-Identifier: Apache-2.0 */
#include <cblas.h>
#include <jemalloc/jemalloc.h>
#include <stdio.h>
#include <time.h>

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

int main() {
	srand(time(NULL));
	/* B, C cols */
	const size_t n = 2000;
	/* A, C rows */
	const size_t m = 2000;
	/* A cols, B rows */
	const size_t k = 3;

	matrix a, b, c;
	matrix_ccreate(&a, n, m, malloc);
	matrix_ccreate(&b, n, m, malloc);
	matrix_ccreate(&c, n, m, malloc);

	generate_test(&a, &b, n);

	printf(";; OpenBLAS res\n");

	clock_t start = clock();
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1,
	            a.data->data, k, b.data->data, n, 0, c.data->data, n);
	clock_t end = clock();
	double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

	// matrix_print(&c);
	printf("Time: %f sec\n", time_spent);
	printf(";;\n\n");

	// matrix_copy

	printf(";; mm res\n");
	start = clock();
	matrix_on_matrix(&a, &b, &c);
	end = clock();

	time_spent = (double)(end - start) / CLOCKS_PER_SEC;
	// matrix_print(&c);
	printf("Time: %f sec\n", time_spent);
	printf(";;\n");

	matrix_destroy(&a, &b, &c);
	return 0;
}
