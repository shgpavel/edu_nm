/* SPDX-License-Identifier: Apache-2.0 */

#include <jemalloc/jemalloc.h>
#include <math.h>
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mult.h"

int generate_test(matrix *a, matrix *b, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			matrix_val(a, i, j) =
			    ((double)rand() / (double)RAND_MAX) - 0.5;
			matrix_val(b, i, j) =
			    ((double)rand() / (double)RAND_MAX) - 0.5;
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
	size_t const alignment = 64;
	if (size % alignment != 0) {
		size = (size + alignment - 1) & ~(alignment - 1);
	}
	return aligned_alloc(alignment, size);
}

void *arealloc(void *data, size_t size) {
	return realloc(data, size);
}

typedef void (*matrix_mult_func)(matrix const *, matrix const *,
                                 matrix const *);
struct stats {
	double avg;
	double p95;
	double p1;
};

/*
struct stats bench(matrix_mult_func multer, size_t n, size_t tests) {
        matrix a, b, c;
        matrix_ccreate(n, n, aalloc, arealloc, &a, &b, &c);

        generate_test(&a, &b, n);

        vector times;
        vector_ccreate(&times, tests, aalloc);

        for (size_t i = 0; i < tests; ++i) {
                clock_t start = clock();
                multer(&a, &b, &c);

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
        matrix_destroy(&a, &b, &c);
        vector_destroy(&times);
        return res;
}
*/
void gen_name(char *buf, size_t size) {
	char const charset[] = "abcdefghijklmnopqrstuvwxyz0123456789";
	for (size_t i = 0; i < size - 1; i++) {
		buf[i] = charset[rand() % (sizeof(charset) - 1)];
	}
	buf[size - 1] = '\0';
}

void test_correctness() {
	size_t n;
	scanf("%zu", &n);

	matrix a, b, c, btr, c2;
	matrix_ccreate(n, n, aalloc, arealloc, &a, &b, &c, &btr, &c2);

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

	matrix_mult_mkl(&a, &btr, &c2);
	matrix_print(&c2);
	printf("\n");

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			int k = fabs(matrix_val(&c, i, j) -
			             matrix_val(&c2, i, j)) < 1e-5;
			if (j != n - 1)
				printf("%d ", k);
			else
				printf("%d", k);
		}
		printf("\n");
	}
	printf("\n");

	matrix_destroy(&a, &b, &c, &btr, &c2);
}

int main() {
	test_correctness();

	/*
	srand(time(NULL));

	char test1[10];
	gen_name(test1, 5);
	char test0[5] = ".csv";
	strcat(test1, test0);

	char test2[15] = "test-";
	strcat(test2, test1);
	char test3[25] = "results/";
	strcat(test3, test2);

	printf(";;\n  results in %s\n;;\n", test3);

	FILE *csvout = fopen(test3, "w");
	if (!csvout) {
	        return -1;
	}

	// TODO activate/de tests otf
	int err = 0;
	err = fprintf(csvout,
	              "n,"
	              "mklavg,mkl95p,mkl1p,"
	              //"naiveavg,naive95p,naive1p,"
	              "resavg,res95p,res1p\n");
	if (err < 0) {
	        return -1;
	}

	size_t const nmin = 64, nmax = 3090;
	size_t const aclen = 18;
	size_t acnum[] = {20, 15, 10, 9, 8, 7, 6, 5, 4,
	                  3,  2,  1,  1, 1, 1, 1, 1, 1};

	for (size_t n = nmin; n < nmax; n += nmin) {
	        double izedn =
	            ((double)n - (double)nmin) / ((double)nmax - (double)nmin);
	        size_t kidx =
	            (size_t)round((double)izedn * ((double)aclen - 1.0));
	        size_t ktests = acnum[kidx];

	        printf("%zu %zu %.3f%%\n", n, ktests, izedn * 100);

	        struct stats resMKL = bench(matrix_mult_mkl, n, ktests);
	        // struct stats resN = bench(matrix_mult_naive, n, ktests);
	        struct stats resR = bench(matrix_mult_fast, n, ktests);

	        fprintf(csvout,
	                "%zu,"
	                "%lg,%lg,%lg,"
	                //"%lg,%lg,%lg,"
	                "%lg,%lg,%lg\n",
	                n, resMKL.avg, resMKL.p95, resMKL.p1,
	                // resN.avg, resN.p95, resN.p1,
	                resR.avg, resR.p95, resR.p1);

	        n += acnum[aclen - 1 - kidx] * nmin;
	}

	err = fclose(csvout);
	if (err) {
	        return -2;
	}
	*/
	return 0;
}
