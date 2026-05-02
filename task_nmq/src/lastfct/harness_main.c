#include <stdio.h>
#include <stdlib.h>

extern int contest_last_nonzero_digit_factorial(int n);

/* QEMU/NMC4 cannot read stdin reliably, so pass inputs via argv:
 *   nmc-qemu solution.abs <n> [iterations]
 */
int main(int argc, char** argv)
{
    int n;
    int i;
    int iters;
    int result;
    volatile int sink;
    char* endptr;

    if (argc < 2) {
        fprintf(stderr, "Usage: nmc-qemu solution.abs <n> [iterations]\n");
        return 1;
    }

    n = (int)strtol(argv[1], &endptr, 10);
    if (*endptr != '\0' || n < 0) {
        fprintf(stderr, "Invalid n: %s\n", argv[1]);
        return 1;
    }

    iters = (argc >= 3) ? atoi(argv[2]) : 1;
    if (iters < 1)
        iters = 1;

    sink = 0;
    for (i = 0; i < iters; i++)
        sink = contest_last_nonzero_digit_factorial(n);

    result = sink;
    printf("%d\n", result);
    return 0;
}
