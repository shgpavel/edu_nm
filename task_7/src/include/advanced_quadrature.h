/* SPDX-License-Identifier: Apache-2.0 */

#ifndef ADVANCED_QUADRATURE
#define ADVANCED_QUADRATURE

#include <stddef.h>

double newton_cotes(double, double);
double gauss_quad(double, double);
double newton_cotes_nxn(size_t, double, double);

#endif
