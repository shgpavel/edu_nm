#ifndef PARTITIONS_H
#define PARTITIONS_H

#include "func.h"
#include "vector.h"

inline vector linspace(size_t num_of_points) {
	vector res;
	vector_init(&res, num_of_points, sizeof(double));
	for (size_t i = 0; i < num_of_points; ++i) {
		double point = A + (double)i * (B - A) / (double) (num_of_points - 1);
		vector_push(&res, &point);
	}
	return res;
}

#endif
