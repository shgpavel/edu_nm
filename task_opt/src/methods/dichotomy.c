/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>
#include <math.h>

#include "func.h"

double dichotomy_method(double a, double b, double epsilon, double delta) {
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}

	double mid = (a + b) / 2;
	double c = mid - delta / 2;
	double d = mid + delta / 2;
	
	
	
	return (a + b) / 2;
}
