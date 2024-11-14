/* SPDX-License-Identifier: Apache-2.0 */

/*
 * Memory Manager Prototype for Numerical Methods (nm)
 *
 * This system consists of two processes
 *
 * 1. Computational Process (C)
 *  - Performs numerical calculations.
 *  - Initially allocates all data on the stack.
 *  - Updates shared memory with current memory usage.
 *  - Upon receiving a signal from the O, switches
 *    to heap allocation to prevent stack overflow.
 *
 * 2. Observer Process (O)
 *  - Monitors C's memory usage via shared memory.
 *  - Sends a signal ( SIGUSR1 ) to C when stack usage
 *    approaches the limit, prompting C to allocate on the heap.
 *
 * NOTE: This thing built 100% for fun, it is only around 20-50% possibility
 *       that you will be able to run this on GNU/Linux, because
 *       stack size behavior is not strict and standardized.
 *       I tested POSIX functions to change stack limits in runtime,
 *       and they are useless on my system, however you can change
 *       stack size by ulimit -s [size] command in the shell.
 *       Remember that stack also contains binary i. e. .code.
 */

#include <fcntl.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <semaphore.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <errno.h>

#define SHM_NAME "/shared_control"

struct control {
	sem_t semaphore;
	uint64_t stack_size;
	uint64_t stack_max;
	int fd;
};

size_t stack_size(char *filename) {
	struct rlimit limit;
	if (getrlimit(RLIMIT_STACK, &limit) != 0) {
		perror("getrlimit");
		return 0;
	}

	FILE *code = fopen(filename, "rb");
	if (code == NULL) {
		perror("open");
		return 0;
	}

	if (fseek(code, 0, SEEK_END) != 0) {
		perror("fseek");
		goto failure;
	}

	long file_size = ftell(code);
	if (file_size == -1) {
		perror("ftell");
		goto failure;
	}

	return (size_t)limit.rlim_cur - (size_t)file_size;
	
 failure:
	fclose(code);
	return 0;
}

struct control *shared_control_create(char *filename) {
	int fd = shm_open(SHM_NAME, O_CREAT | O_RDWR, 0666);
	if (fd == -1) {
		perror("shm_open");
		return NULL;
	}

	if (ftruncate(fd, sizeof(struct control)) == -1) {
		perror("ftruncate");
		goto failure;
	}

	struct control *shared_data = mmap(NULL, sizeof(struct control),
																		 PROT_READ | PROT_WRITE,
																		 MAP_SHARED, fd, 0);

	if (shared_data == MAP_FAILED) {
		perror("mmap");
		goto failure;
	}

	if (sem_init(&shared_data->semaphore, 1, 1) == -1) {
		perror("sem_init");
		munmap(shared_data, sizeof(struct control));
		goto failure;
	}

	shared_data->stack_size = 0;
	shared_data->stack_max = stack_size(filename);
	shared_data->fd = fd;
	return shared_data;
	
 failure:
	close(fd);
	return NULL;
}

struct control *shared_control_open() {
	int fd = shm_open(SHM_NAME, O_RDWR, 0666);
	if (fd == -1) {
		perror("shm_open");
		return NULL;
	}

	struct control *shared_data = mmap(NULL, sizeof(struct control),
																		 PROT_READ | PROT_WRITE,
																		 MAP_SHARED, fd, 0);
	if (shared_data == MAP_FAILED) {
		perror("mmap");
		close(fd);
		return NULL;
	}

	if (shared_data->fd) {
		
	} else {
		close(fd);
	}
	return shared_data; 
}

void shared_control_destroy(struct control *shared_data, int *shm_fd) {
	if (munmap(shared_data, sizeof(struct control)) == -1) {
		perror("munmap");
		exit(EXIT_FAILURE);
	}

	if (close(*shm_fd) == -1) {
		perror("close");
		exit(EXIT_FAILURE);
	}
}

void o() {
	for (size_t i = 0; i < 1000; ++i) {
		printf("1\t");
	}
}

void c() {
	for (size_t i = 0; i < 1000; ++i) {
		printf("2\t");
	}
}

int main(int argc, char *argv[]) {
	
	pid_t pid = fork();

	if (pid == 0) {
		printf("Child process: PID=%d\n", getpid());
		c();
	} else if (pid > 0) {
		printf("Parent process: PID=%d\n", getpid());
		o();
	} else {
		printf("Error: Could not create child\n");
		return -1;
	}
	printf("\n");
	
	return 0;
}
