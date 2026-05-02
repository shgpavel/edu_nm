#include <stdio.h>
#include <stdlib.h>

typedef unsigned int u32;

extern u32 contest_pow_mod(u32 base, u32 exp, u32 mod);

/* QEMU/NMC4 cannot read stdin reliably, so pass inputs via argv:
 *   nmc-qemu solution.abs <base> <exp> <mod> [iterations]
 */
int main(int argc, char** argv)
{
    int i;
    int iters;
    u32 base;
    u32 exp;
    u32 mod;
    u32 result;
    volatile u32 sink;
    char* endptr;

    if (argc < 4) {
        fprintf(stderr, "Usage: nmc-qemu solution.abs <base> <exp> <mod> [iterations]\n");
        return 1;
    }

    base = (u32)strtoul(argv[1], &endptr, 0);
    if (*endptr != '\0') {
        fprintf(stderr, "Invalid base: %s\n", argv[1]);
        return 1;
    }

    exp = (u32)strtoul(argv[2], &endptr, 0);
    if (*endptr != '\0') {
        fprintf(stderr, "Invalid exp: %s\n", argv[2]);
        return 1;
    }

    mod = (u32)strtoul(argv[3], &endptr, 0);
    if (*endptr != '\0') {
        fprintf(stderr, "Invalid mod: %s\n", argv[3]);
        return 1;
    }

    iters = (argc >= 5) ? atoi(argv[4]) : 1;
    if (iters < 1)
        iters = 1;

    sink = 0U;
    for (i = 0; i < iters; i++)
        sink = contest_pow_mod(base, exp, mod);

    result = sink;
    printf("%u\n", result);
    return 0;
}
