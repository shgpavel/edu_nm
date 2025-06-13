/* SPDX-License-Identifier: Apache-2.0 */

#ifndef COMPOUND_H
#define COMPOUND_H

#include <stddef.h>

typedef double (*quad_func)(double, double);
typedef double (*quad_func_any_n)(size_t, double, double);

double compound(size_t, quad_func);
double compound_ncwgau(size_t, size_t, quad_func_any_n, quad_func);

#endif
