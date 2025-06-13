/* SPDX-License-Identifier: Apache-2.0 */

#include "richardson.h"

#include <jemalloc/jemalloc.h>
#include <math.h>
#include <stdio.h>

#include "aitken.h"
#include "compound.h"
#include "func.h"
#include "matrix.h"
#include "qr_avx.h"
#include "vector.h"

struct r_data richardson_step(int m, vector *sums, int n, double first_step) {
	struct r_data result = {nan(""), nan(""), nan(""), nan("")};

	matrix t;
	matrix_cinit(&t, n, n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		matrix_val(&t, i, 0) = 1.0;
		double base = first_step * pow(0.5, i);
		for (int j = 0; j < n - 1; ++j) {
			matrix_val(&t, i, j + 1) = pow(base, j + m + 1);
		}
	}

	vector b;
	vector_init_copy(&b, sums);

	vector *v = qr_avx(&t, &b);

	result.value = vector_val(v, 0);
	result.error = fabs(vector_val(sums, n - 1) - result.value);
	/* last elem */
	result.m = vector_val(v, n - 1);

	vector_free(v);
	free(v);

	vector_free(&b);
	matrix_free(&t);

	return result;
}

struct r_data richardson(quad_func qf, double eps, size_t initial_steps) {
	struct r_data result = {nan(""), nan(""), nan(""), nan("")};

	vector sums;
	vector_cinit(&sums, 2, sizeof(double));

	vector_val(&sums, 0) = compound(initial_steps, qf);
	vector_val(&sums, 1) = compound(2 * initial_steps, qf);

	struct r_data *report =
	    (struct r_data *)malloc(10 * sizeof(struct r_data));

	double first_step = B - A / initial_steps;

	report[0].step = first_step;
	report[0].m = 0;
	report[0].value = vector_val(&sums, 0);
	report[0].error = 0;

	report[1].step = (B - A) / (2 * initial_steps);
	report[1].m = 0;
	report[1].value = vector_val(&sums, 1);
	report[1].error = 0;

	// printf("sums init: %zu\n", sums.size);
	// vector_print(&sums);

	int max_iterations = 100, report_len = 2, report_cap = 10;
	for (int iter = 0; iter < max_iterations; ++iter) {
		/* next_n = 2^(sums_len) * initial_steps
		         bc 2^(sums_len) == (1 << sums_len) */
		int next_n = initial_steps * (1 << sums.size);
		double new_sum = compound(next_n, qf);

		vector_push(&sums, &new_sum);

		if (sums.size < 3) {
			continue;
		}
		// printf("ss: %zu\n", sums.size);
		// vector_print(&sums);

		double last_three[3];
		last_three[0] = vector_val(&sums, sums.size - 3);
		last_three[1] = vector_val(&sums, sums.size - 2);
		last_three[2] = vector_val(&sums, sums.size - 1);
		double m_val = aitken_process(last_three, 2);
		int m_int = (int)m_val;

		struct r_data rich_res =
		    richardson_step(m_int, &sums, sums.size, first_step);
		if (isnan(rich_res.value)) {
			break;
		}

		if (report_len >= report_cap) {
			report_cap *= 2;
			report = (struct r_data *)realloc(
			    report, report_cap * sizeof(struct r_data));
		}
		report[report_len].step = first_step * pow(0.5, sums.size);
		report[report_len].m = m_val;
		report[report_len].value = rich_res.value;
		report[report_len].error = rich_res.error;
		report_len++;

		if (rich_res.error < eps) {
			printf("Richardson success\n");
			result.value = rich_res.value;
			result.m = m_val;
			result.error = rich_res.error;
			result.step = report[report_len].step;
			goto clear;
		}
	}
	printf("Richardson failed\n");

clear:
	free(report);
	vector_free(&sums);
	return result;
}

struct r_data step_guessing(quad_func qf, double eps) {
	struct r_data result;
	vector sums;
	vector_cinit(&sums, 3, sizeof(double));

	vector_val(&sums, 0) = compound(1, qf);
	vector_val(&sums, 1) = compound(2, qf);
	vector_val(&sums, 2) = compound(4, qf);

	vector reports;
	vector_init(&reports, 2, sizeof(struct r_data));

	struct r_data first_rep = {.error = nan(""),
	                           .step = B - A,
	                           .m = nan(""),
	                           .value = vector_val(&sums, 0)};
	struct r_data sec_rep = {.error = nan(""),
	                         .step = (B - A) / 2,
	                         .m = nan(""),
	                         .value = vector_val(&sums, 1)};

	vector_push(&reports, &first_rep);
	vector_push(&reports, &sec_rep);

	double aip[3];
	aip[0] = vector_val(&sums, sums.size - 3);
	aip[1] = vector_val(&sums, sums.size - 2);
	aip[2] = vector_val(&sums, sums.size - 1);

	double m = aitken_process(aip, 2);

	struct r_data thr_rep =
	    richardson_step((int)m, &sums, sums.size, B - A);
	struct r_data fr_rep = {.error = thr_rep.error,
	                        .step = (B - A) / 4,
	                        .m = thr_rep.m,
	                        .value = vector_val(&sums, 4)};

	if (fr_rep.error < eps) {
		result = fr_rep;
		goto clear;
	}

	double base = eps / thr_rep.m;
	double approx_h = pow(base, 1.0 / m);
	size_t initial_steps = (B - A) / approx_h;
	if (initial_steps > 100) {
		initial_steps = 5;
	}
	// printf("%zu\n", initial_steps);

	result = richardson(qf, eps, initial_steps);

clear:
	vector_free(&sums);
	vector_free(&reports);

	return result;
}
