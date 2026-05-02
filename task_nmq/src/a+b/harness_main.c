#include <stdio.h>
#include <stdlib.h>
#include <time.h>

extern int contest_a_plus_b(int a, int b);

/* NMC4 has CHAR_BIT=32, so stdin/file scanf is broken under QEMU semihosting.
   Pass inputs as argv: nmc-qemu solution.abs <a> <b> [iterations]
   Default iterations: 1000000 */
int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: nmc-qemu solution.abs <a> <b> [iterations]\n");
        return 1;
    }

    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    int iters = (argc >= 4) ? atoi(argv[3]) : 1000000;

    /* Warmup */
    volatile int sink = contest_a_plus_b(a, b);

    clock_t start = clock();
    int i;
    for (i = 0; i < iters; i++) {
        sink = contest_a_plus_b(a, b);
    }
    clock_t end = clock();

    long long ticks = (long long)(end - start);
    long long cps   = (long long)CLOCKS_PER_SEC;

    printf("%d\n", sink);
    fprintf(stderr, "result:       %d\n",   contest_a_plus_b(a, b));
    fprintf(stderr, "iterations:   %d\n",   iters);
    fprintf(stderr, "total_ticks:  %lld\n", ticks);
    fprintf(stderr, "clock_per_sec:%lld\n", cps);
    if (ticks > 0)
        fprintf(stderr, "ns_per_call:  %lld\n", ticks * 1000000000LL / cps / iters);
    else
        fprintf(stderr, "ns_per_call:  <1 tick (increase iterations)\n");
    return 0;
}
