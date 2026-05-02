#include <stdio.h>
#include <stdlib.h>

#include "contest_matmul_i32.h"

/*
 * QEMU harness — stdin is broken on NMC4/QEMU, pass inputs via argv.
 *
 * Mode 1 — explicit matrix (correctness / open tests):
 *   solution.abs <m> <k> <n> <a[0]> ... <a[m*k-1]> <b[0]> ... <b[k*n-1]>
 *
 * Mode 2 — generated matrix benchmark (no print, just timing):
 *   solution.abs gen <m> <k> <n> <iters>
 *   Fills A and B with deterministic pseudo-random values, runs iters times.
 *   Output: nothing (time from host shell).
 *
 * Mode 3 — generated correctness self-check:
 *   solution.abs checkgen <m> <k> <n>
 *   Fills A and B deterministically, runs contest_matmul_i32 once, then
 *   compares against a local scalar reference and returns nonzero on mismatch.
 */

#define MAX_DIM  1024
#define MAX_ELEMS (MAX_DIM * MAX_DIM)

/* Place large arrays in EMI (external memory) — NMMB1 (.bss) is only 64K words */
static int A[MAX_ELEMS] __attribute__((section(".bss.emi")));
static int B[MAX_ELEMS] __attribute__((section(".bss.emi")));
static int C[MAX_ELEMS] __attribute__((section(".bss.emi")));

/* Simple LCG for deterministic fill without stdlib rand (which may call Mul32) */
static int lcg(int* state) {
    *state = (*state) * 1664525 + 1013904223;
    return (*state) >> 16;   /* use upper 16 bits */
}

static int fill_value(int* state) {
    return (lcg(state) & 7) - 3;   /* bounded values keep reference sums small */
}

static int check_generated_case(int m, int k, int n)
{
    int seed = 42;
    int i, row, col, t;

    for (i = 0; i < m * k; i++) A[i] = fill_value(&seed);
    for (i = 0; i < k * n; i++) B[i] = fill_value(&seed);

    contest_matmul_i32(A, B, C, m, k, n);

    for (row = 0; row < m; row++) {
        for (col = 0; col < n; col++) {
            int expected = 0;
            for (t = 0; t < k; t++)
                expected += A[row * k + t] * B[t * n + col];
            if (C[row * n + col] != expected) {
                fprintf(stderr,
                        "Mismatch at (%d,%d): expected %d got %d\n",
                        row, col, expected, C[row * n + col]);
                return 2;
            }
        }
    }

    return 0;
}

int main(int argc, char** argv)
{
    /* ---- Mode 3: generated correctness self-check ---- */
    if (argc >= 2 &&
        argv[1][0] == 'c' && argv[1][1] == 'h' && argv[1][2] == 'e' &&
        argv[1][3] == 'c' && argv[1][4] == 'k' && argv[1][5] == 'g' &&
        argv[1][6] == 'e' && argv[1][7] == 'n' && argv[1][8] == '\0') {
        if (argc < 5) {
            fprintf(stderr, "Usage: solution.abs checkgen m k n\n");
            return 1;
        }
        return check_generated_case(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    }

    /* ---- Mode 4: multibench — one QEMU invocation for many shapes ---- */
    /*    solution.abs multi iters m1 k1 n1 m2 k2 n2 ...                    */
    if (argc >= 6 && argv[1][0] == 'm' && argv[1][1] == 'u' && argv[1][2] == 'l' &&
        argv[1][3] == 't' && argv[1][4] == 'i' && argv[1][5] == '\0') {
        int iters = atoi(argv[2]);
        int shape_argc = argc - 3;
        int nshapes = shape_argc / 3;
        int s, i, seed = 42, fill_len;
        /* Fill A, B once with max extent needed — use first shape's K*max N */
        for (i = 0; i < MAX_ELEMS; i++) A[i] = lcg(&seed);
        for (i = 0; i < MAX_ELEMS; i++) { B[i] = lcg(&seed); if (i > 65535) break; }
        (void)fill_len;
        for (s = 0; s < nshapes; s++) {
            int m = atoi(argv[3 + 3*s]);
            int k = atoi(argv[3 + 3*s + 1]);
            int n = atoi(argv[3 + 3*s + 2]);
            for (i = 0; i < iters; i++)
                contest_matmul_i32(A, B, C, m, k, n);
        }
        return 0;
    }

    /* ---- Mode 2: generated benchmark ---- */
    if (argc >= 2 && argv[1][0] == 'g' && argv[1][1] == 'e' && argv[1][2] == 'n') {
        if (argc < 6) {
            fprintf(stderr, "Usage: solution.abs gen m k n iters\n");
            return 1;
        }
        int m     = atoi(argv[2]);
        int k     = atoi(argv[3]);
        int n     = atoi(argv[4]);
        int iters = atoi(argv[5]);

        int seed = 42, i;
        for (i = 0; i < m * k; i++) A[i] = lcg(&seed);
        for (i = 0; i < k * n; i++) B[i] = lcg(&seed);

        for (i = 0; i < iters; i++)
            contest_matmul_i32(A, B, C, m, k, n);
        return 0;
    }

    /* ---- Mode 1: explicit matrix ---- */
    if (argc < 4) {
        fprintf(stderr, "Usage: solution.abs m k n a00 ... b00 ...\n"
                        "       solution.abs checkgen m k n\n"
                        "       solution.abs gen m k n iters\n");
        return 1;
    }

    int m = atoi(argv[1]);
    int k = atoi(argv[2]);
    int n = atoi(argv[3]);

    int need = 4 + m * k + k * n;
    if (argc < need) {
        fprintf(stderr, "Need %d args, got %d\n", need, argc);
        return 1;
    }

    int i;
    for (i = 0; i < m * k; i++)
        A[i] = atoi(argv[4 + i]);
    for (i = 0; i < k * n; i++)
        B[i] = atoi(argv[4 + m * k + i]);

    /* Optional trailing iteration count (suppress output) */
    int iters = 1, print_result = 1;
    if (argc > need) {
        iters = atoi(argv[need]);
        print_result = 0;
    }

    for (i = 0; i < iters; i++)
        contest_matmul_i32(A, B, C, m, k, n);

    if (print_result) {
        int row, col;
        for (row = 0; row < m; row++) {
            for (col = 0; col < n; col++) {
                if (col) printf(" ");
                printf("%d", C[row * n + col]);
            }
            printf("\n");
        }
    }

    return 0;
}
