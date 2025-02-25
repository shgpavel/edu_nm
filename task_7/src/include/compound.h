/* SPDX-License-Identifier: Apache-2.0 */

#ifndef COMPOSITING_H
#define COMPOSITING_H

typedef double (*quad_func)(double, double);
typedef double (*quad_func_any_n)(size_t, double, double);

double compound(size_t, quad_func);
double compound_ncwgau(size_t, size_t, quad_func_any_n, quad_func);

#endif
