/* SPDX-License-Identifier: Apache-2.0 */

#include <math.h>

#define SOME_TOL 1e-6

double aitken_process(double s[3], double l) {
	double upper_val = s[2] - s[1];
	double lower_val = s[1] - s[0];
	if (fabs(lower_val) < SOME_TOL) {
		return nan("");
	}

	double upper_log = upper_val / lower_val;
	if (upper_log < 0) {
		upper_log = log(fabs(upper_log));
	} else {
		upper_log = log(upper_log);
	}
	return -upper_log / log(l);
}
