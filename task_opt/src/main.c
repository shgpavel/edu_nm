/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>

#include "dichotomy.h"
#include "direct_search.h"
#include "func.h"
#include "golden.h"

int main() {
	double const a = 0;
	double const b = 1;

	double const epsilon = 1e-6;
	double const delta = 1e-7;
	double const step = -1e-1;

	char res[2][25];

	printf(";; Golden ratio minimizer\n");
	double x = golden(a, b, epsilon);

	format_hex(res[0], 25, x);
	format_hex(res[1], 25, func(x));

	printf("argmin = %s f(argmin) = %s\n", res[0], res[1]);
	printf(";;\n\n");

	printf(";; Dichotomy minimizer\n");
	x = dichotomy_method(a, b, epsilon, delta);

	format_hex(res[0], 25, x);
	format_hex(res[1], 25, func(x));

	printf("argmin = %s f(argmin) = %s\n", res[0], res[1]);
	printf(";;\n\n");

	printf(";; Direct search minimizer\n");
	x = direct_search(b, step);

	format_hex(res[0], 25, x);
	format_hex(res[1], 25, func(x));

	printf("argmin = %s f(argmin) = %s\n", res[0], res[1]);
	printf(";;\n");
}
