/* SPDX-License-Identifier: Apache-2.0 */

#ifndef PARTITIONS_H
#define PARTITIONS_H

#include "func.h"
#include "vector.h"

inline vector linspace(double a, double b, size_t num_of_points) {
	vector res;
	vector_init(&res, num_of_points, sizeof(double));
	for (size_t i = 0; i < num_of_points; ++i) {
		double point =
		    a + (double)i * (b - a) / (double)(num_of_points - 1);
		vector_push(&res, &point);
	}
	return res;
}

#endif
