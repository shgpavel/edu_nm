/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>

#include "func.h"

int main() {
	/*
	long nprocs = sysconf(_SC_NPROCESSORS_ONLN);
	if (nprocs < 1) {
		perror("sysconf");
		return -1;
	}

	pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * nprocs);
	if (threads == NULL) {
		perror("malloc");
		return -1;
	}

	for (long i = 0; i < nprocs; ++i) {
		if (pthread_create(&threads[i], NULL, thread_function, (void*)i) != 0) {}
	}
	return 0;
	*/
}
