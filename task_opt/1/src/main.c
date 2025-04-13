/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>

#include "dichotomy.h"
#include "direct_search.h"
#include "func.h"
#include "golden.h"

int main() {
	double const a = 0;
	double const b = 1;
	double const epsilon = 1e-4;

	printf(";; Golden ratio minimizer\n");
	double x = golden(a, b, epsilon);
	//printf("x = %.3a f(x) = %.3a\n", x, func(x));
	printf("argmin = %lg f(argmin) = %lg span = %d\n", x, func(x),  0);
	printf(";;\n\n");
	
	printf(";; Dichotomy minimizer\n");
	x = dichotomy_method(a, b, epsilon);
	//printf("x = %a f(x) = %a\n", x, func(x));
	printf("argmin = %lg f(argmin) = %lg span = %d\n", x, func(x), 0);
	printf(";;\n\n");

	printf(";; Direct search minimizer\n");
	x = direct_search(a, b);
	//printf("x = %a f(x) = %a\n", x, func(x));
	printf("argmin = %lg f(argmin) = %lg span = %d\n", x, func(x), 0);
	printf(";;\n");
}
