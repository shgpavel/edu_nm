#include <stdio.h>
#include <stdlib.h>

struct rc4_ctx {
    unsigned char S[256];
    int i;
    int j;
};

extern void contest_rc4_init(struct rc4_ctx *ctx, const unsigned char *key, int key_len);
extern void contest_rc4_crypt(struct rc4_ctx *ctx, const unsigned char *input,
                               unsigned char *output, int len);

static unsigned char key_buf[256];
static unsigned char input_buf[4096];
static unsigned char output_buf[4096];

/*
 * Mode 1 — explicit (correctness / open tests):
 *   solution.abs <key_len> <k0..kN-1> <data_len> <d0..dM-1> [iters]
 *   iters > 1: suppress output (bench mode)
 *
 * Mode 2 — generated benchmark (avoids huge argv):
 *   solution.abs gen <key_len> <data_len> <iters>
 *   Fills key and input with deterministic values (i & 0xFF).
 *   Prints nothing; time from host shell.
 */
int main(int argc, char **argv)
{
    int i;

    if (argc < 2) {
        fprintf(stderr, "Usage: solution.abs key_len k0..kN data_len d0..dM [iters]\n"
                        "       solution.abs gen key_len data_len iters\n");
        return 1;
    }

    /* ---- Mode 2: gen ---- */
    if (argv[1][0] == 'g' && argv[1][1] == 'e' && argv[1][2] == 'n') {
        if (argc < 5) {
            fprintf(stderr, "Usage: solution.abs gen key_len data_len iters\n");
            return 1;
        }
        int key_len  = atoi(argv[2]);
        int data_len = atoi(argv[3]);
        int iters    = atoi(argv[4]);

        for (i = 0; i < key_len;  i++) key_buf[i]   = (unsigned char)(i & 0xFF);
        for (i = 0; i < data_len; i++) input_buf[i]  = (unsigned char)(i & 0xFF);

        struct rc4_ctx ctx;
        for (i = 0; i < iters; i++) {
            contest_rc4_init(&ctx, key_buf, key_len);
            contest_rc4_crypt(&ctx, input_buf, output_buf, data_len);
        }
        return 0;
    }

    /* ---- Mode 1: explicit ---- */
    int idx = 1;
    int key_len  = atoi(argv[idx++]);
    for (i = 0; i < key_len; i++)
        key_buf[i] = (unsigned char)atoi(argv[idx++]);

    int data_len = atoi(argv[idx++]);
    for (i = 0; i < data_len; i++)
        input_buf[i] = (unsigned char)atoi(argv[idx++]);

    int iters = (idx < argc) ? atoi(argv[idx]) : 1;

    struct rc4_ctx ctx;
    for (i = 0; i < iters; i++) {
        contest_rc4_init(&ctx, key_buf, key_len);
        contest_rc4_crypt(&ctx, input_buf, output_buf, data_len);
    }

    if (iters <= 1 || idx >= argc) {
        for (i = 0; i < data_len; i++)
            printf("%d ", output_buf[i]);
        printf("\n");
    }

    return 0;
}
