/* SPDX-License-Identifier: Apache-2.0 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include "advanced_quadrature.h"
#include "aitken.h"
#include "compound.h"
#include "func.h"
#include "richardson.h"
#include "simple_quadrature.h"

int get_data_tonethree() {
	printf(";; main 1.3\n");
	double const val = 6.2389826790;
	double const val2 = 27.5664955364;

	int err = 0;

	// if (mkdir("out", 0755) != 0)
	//   return -1;

	FILE *csvout = fopen("out/tonethree.csv", "w");
	if (!csvout) {
		return -1;
	}

	err = fprintf(csvout,
	              "n, D_left(n), D_mid(n), D_trapz(n), D_simpson(n), "
	              "D_nc(n), D_gs(n), D_ts(n)\n");
	if (err < 0) {
		return -2;
	}

	for (size_t i = 1; i <= 100; ++i) {
		double absdiff[6];
		absdiff[0] = fabs(val - compound(i, quad_left));
		absdiff[1] = fabs(val - compound(i, quad_mid));
		absdiff[2] = fabs(val - compound(i, quad_trapz));
		absdiff[3] = fabs(val - compound(i, quad_simpson));

		absdiff[4] = fabs(val2 - compound(i, newton_cotes));
		absdiff[5] = fabs(val2 - compound(i, gauss_quad));

		err = fprintf(
		    csvout,
		    "%zu, %.15lf, %.15lf, %.15lf, %.15lf, %.15lf, %.15lf\n", i,
		    absdiff[0], absdiff[1], absdiff[2], absdiff[3], absdiff[4],
		    absdiff[5]);
		if (err < 0) {
			return -3;
		}
	}

	err = fclose(csvout);
	if (err) {
		return -4;
	}
	printf("Info: 1.3 saved to `out/ottonethree.csv`\n");
	printf(";;\n\n");
	return 0;
}

int get_opt_step_ttwoone() {
	printf(";; main 2.1 Newton-Cotes Richardson\n");
	struct r_data res = richardson(newton_cotes, 1e-6, 10);
	printf("S = %.7lf e = %e m = %.7lf on step = %e\n", res.value,
	       res.error, res.m, res.step);
	printf(";;\n\n");

	printf(";; main 2.2 Gauss Richardson\n");
	struct r_data res2 = richardson(gauss_quad, 1e-6, 3);
	printf("S = %.7lf e = %e m = %.7lf on step = %e\n", res2.value,
	       res2.error, res2.m, res2.step);
	printf(";;\n\n");
	return 0;
}

int guessing_ttwothree() {
	printf(";; main 2.3 Newton-Cotes Guessing\n");
	struct r_data res = step_guessing(newton_cotes, 1e-6);
	printf("S = %.7lf e = %e m = %.7lf on step = %e\n", res.value,
	       res.error, res.m, res.step);
	printf(";;\n\n");
	return 0;
}

int check_newton_nxn_penone() {
	printf(";; pen 1.1 Newton-Cotes nxn\n");
	printf("i\tS\n");
	for (size_t i = 2; i < 8; ++i) {
		printf("%zu\t%.9lf\n", i, newton_cotes_nxn(i, A, B));
	}
	printf(";;\n\n");
	return 0;
}

int get_data_pentwo() {
	const double val2 = 27.5664955364;
	printf(";; pen.2\n");
	int err = 0;

	FILE *csvout = fopen("out/pentwo.csv", "w");
	if (!csvout) {
		return -1;
	}

	err = fprintf(csvout, "n, D_nc_compound(n), D_nc_def(2n + 1)\n");
	if (err < 0) {
		return -2;
	}

	for (size_t i = 1; i <= 30; ++i) {
		double absdiff_ncc = fabs(val2 - compound(i, newton_cotes));
		double absdiff_ncg =
		    fabs(val2 - newton_cotes_nxn(2 * i + 1, A, B));

		err = fprintf(csvout, "%zu, %.15lf, %.15lf\n", i, absdiff_ncc,
		              absdiff_ncg);
		if (err < 0) {
			return -3;
		}
	}

	err = fclose(csvout);
	if (err) {
		return -4;
	}

	printf("Info: pen.2 saved to `out/pentwo.csv`\n");
	printf(";;\n\n");
	return 0;
}

int compound_extra_penthree() {
	printf(";; pen.34 Newton-Cotes3 + Gauss3\n");
	printf("L\tS\tM\n");

	double s[3] = {0, 0, 0};
	for (size_t l = 1; l <= 50; ++l) {
		double sval =
		    compound_ncwgau(l, 3, newton_cotes_nxn, gauss_quad);

		printf("%zu\t%.9lf\t", l, sval);
		s[0] = s[1];
		s[1] = s[2];
		s[2] = sval;

		if (l >= 3) {
			printf("%.6lf\n", aitken_process(s, l));
		} else {
			printf("%lf\n", nan(""));
		}
	}
	printf(";;\n\n");

	printf(";; pen.34 Newton-Cotes6 + Gauss3\n");
	printf("L\tS\tM\n");

	for (size_t l = 1; l <= 50; ++l) {
		double sval =
		    compound_ncwgau(l, 6, newton_cotes_nxn, gauss_quad);
		// int cond = (floor(log2(l)) == log2(l));

		printf("%zu\t%.9lf\t", l, sval);

		s[0] = s[1];
		s[1] = s[2];
		s[2] = sval;

		if (l >= 3) {
			printf("%.6lf\n", aitken_process(s, l));
		} else {
			printf("%lf\n", nan(""));
		}
	}
	printf(";;\n\n");
	return 0;
}

int look_basic_funcs_penfive() {
	const double eps = 1e-4;

	printf(";; pen.5 Left-Rect\n");
	struct r_data res = richardson(quad_left, eps, 1);
	printf("S = %.7lf e = %e m = %.7lf on step = %e\n", res.value,
	       res.error, res.m, res.step);
	printf(";;\n\n");

	printf(";; pen.5 Mid-Rect\n");
	struct r_data res2 = richardson(quad_left, eps, 1);
	printf("S = %.7lf e = %e m = %.7lf on step = %e\n", res2.value,
	       res2.error, res2.m, res2.step);
	printf(";;\n\n");

	printf(";; pen.5 Trapz\n");
	struct r_data res3 = richardson(quad_trapz, eps, 1);
	printf("S = %.7lf e = %e m = %.7lf on step = %e\n", res3.value,
	       res3.error, res3.m, res3.step);
	printf(";;\n\n");

	printf(";; pen.5 Simpson\n");
	struct r_data res4 = richardson(quad_simpson, eps, 1);
	printf("S = %.7lf e = %e m = %.7lf on step = %e\n", res4.value,
	       res4.error, res4.m, res4.step);
	printf(";;\n\n");
	return 0;
}

int main() {
	srand(time(NULL));

	int err = 0;
	/* 1.3 */
	err = get_data_tonethree();
	if (err != 0) {
		fprintf(stderr, "Error: 1.3\n");
		return -1;
	}
	/* 2.1 2.2 */
	get_opt_step_ttwoone();

	/* 2.3 */
	guessing_ttwothree();

	/* 1 */
	err = check_newton_nxn_penone();

	/* 2 */
	err = get_data_pentwo();

	/* 34 */
	err = compound_extra_penthree();

	/* 5 */
	err = look_basic_funcs_penfive();

	return 0;
}
