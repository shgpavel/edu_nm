/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>

#include "golden.h"
#include "func.h"
#include "include/dichotomy.h"

int main() {
	const double a = 0;
	const double b = 1;
	const double epsilon = 1e-4;
	const double delta = 1e-2;
	
	printf(";; Golden ratio minimizer\n");
	double x = golden(a, b, epsilon);
	printf("x = %.3a f(x) = %.3a itr = %zu\n", x, func(x), 1);
	printf(";;\n\n");

	printf(";; Dichotomy minimizer\n");
	x = dichotomy_method(a, b, epsilon, delta);
	printf("x = %a f(x) = %a itr = %zu\n", x, func(x), 1);
	printf(";;\n\n");
}
