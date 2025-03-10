/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>

#include "dichotomy.h"
#include "direct_search.h"
#include "func.h"
#include "golden.h"

int main() {
	const double a = 0;
	const double b = 1;
	const double epsilon = 1e-4;
	const double delta = 1e-2;

	printf(";; Golden ratio minimizer\n");
	double x = golden(a, b, epsilon);
	printf("x = %.3a f(x) = %.3a\n", x, func(x));
	printf(";;\n\n");

	/*
	printf(";; Dichotomy minimizer\n");
	x = dichotomy_method(a, b, epsilon, delta);
	printf("x = %a f(x) = %a\n", x, func(x));
	printf(";;\n\n");
	*/

	printf(";; Direct search minimizer\n");
	x = direct_search(a, b, epsilon);
	printf("x = %a f(x) = %a\n", x, func(x));
	printf(";;\n");
}
