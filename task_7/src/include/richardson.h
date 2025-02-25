/* SPDX-License-Identifier: Apache-2.0 */

#ifndef RICHARDSON_H
#define RICHARDSON_H

#include <stddef.h>

#include "compound.h"
#include "vector.h"

struct r_data {
	double step;
	double m;
	double value;
	double error;
};

struct r_data richardson(quad_func, double, size_t);
struct r_data richardson_step(int, vector *, int, double);

struct r_data step_guessing(quad_func, double);

#endif
