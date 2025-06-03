/* SPDX-License-Identifier: Apache-2.0 */

#ifndef FUNC_H
#define FUNC_H

#include <math.h>
#include <string.h>

#define VDIGIT 6

static inline double func(double x) {
	return exp(-2.6 * x) + exp(1.7 * x) + 0.6 * x;
}

static void format_hex(char *out, size_t len, double val) {
	char tmp[25];
	snprintf(tmp, sizeof(tmp), "%.*a", VDIGIT, val);
	char *p = strstr(tmp, "0x");
	if (p) memmove(p, p + 2, strlen(p + 2) + 1);

	strncpy(out, tmp, len);

	if (len > 0) out[len - 1] = '\0';
}

#endif
