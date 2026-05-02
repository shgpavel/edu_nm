#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define FFT_N 128
#define FFT_WORDS (FFT_N * FFT_N)
#define FFT_DC_SCALE (FFT_N * FFT_N)

typedef struct {
    float re;
    float im;
} complex_float;

void contest_fft_2d(complex_float *matrix, const float *twiddle_cos, const float *twiddle_sin);

static const float k_base_cos[64] = {
    1.000000000f, 0.998795456f, 0.995184727f, 0.989176510f,
    0.980785280f, 0.970031253f, 0.956940336f, 0.941544065f,
    0.923879533f, 0.903989293f, 0.881921264f, 0.857728610f,
    0.831469612f, 0.803207531f, 0.773010453f, 0.740951125f,
    0.707106781f, 0.671558955f, 0.634393284f, 0.595699304f,
    0.555570233f, 0.514102744f, 0.471396737f, 0.427555093f,
    0.382683432f, 0.336889853f, 0.290284677f, 0.242980180f,
    0.195090322f, 0.146730474f, 0.098017140f, 0.049067674f,
    0.000000000f, -0.049067674f, -0.098017140f, -0.146730474f,
    -0.195090322f, -0.242980180f, -0.290284677f, -0.336889853f,
    -0.382683432f, -0.427555093f, -0.471396737f, -0.514102744f,
    -0.555570233f, -0.595699304f, -0.634393284f, -0.671558955f,
    -0.707106781f, -0.740951125f, -0.773010453f, -0.803207531f,
    -0.831469612f, -0.857728610f, -0.881921264f, -0.903989293f,
    -0.923879533f, -0.941544065f, -0.956940336f, -0.970031253f,
    -0.980785280f, -0.989176510f, -0.995184727f, -0.998795456f,
};

static const float k_base_sin[64] = {
    0.000000000f, 0.049067674f, 0.098017140f, 0.146730474f,
    0.195090322f, 0.242980180f, 0.290284677f, 0.336889853f,
    0.382683432f, 0.427555093f, 0.471396737f, 0.514102744f,
    0.555570233f, 0.595699304f, 0.634393284f, 0.671558955f,
    0.707106781f, 0.740951125f, 0.773010453f, 0.803207531f,
    0.831469612f, 0.857728610f, 0.881921264f, 0.903989293f,
    0.923879533f, 0.941544065f, 0.956940336f, 0.970031253f,
    0.980785280f, 0.989176510f, 0.995184727f, 0.998795456f,
    1.000000000f, 0.998795456f, 0.995184727f, 0.989176510f,
    0.980785280f, 0.970031253f, 0.956940336f, 0.941544065f,
    0.923879533f, 0.903989293f, 0.881921264f, 0.857728610f,
    0.831469612f, 0.803207531f, 0.773010453f, 0.740951125f,
    0.707106781f, 0.671558955f, 0.634393284f, 0.595699304f,
    0.555570233f, 0.514102744f, 0.471396737f, 0.427555093f,
    0.382683432f, 0.336889853f, 0.290284677f, 0.242980180f,
    0.195090322f, 0.146730474f, 0.098017140f, 0.049067674f,
};

static complex_float g_matrix[FFT_WORDS] __attribute__((section(".bss.emi")));
static complex_float g_template[FFT_WORDS] __attribute__((section(".bss.emi")));
static complex_float g_expected[FFT_WORDS] __attribute__((section(".bss.emi")));
static complex_float g_ref_scratch[FFT_WORDS] __attribute__((section(".bss.emi")));
static complex_float g_ref_column[FFT_N] __attribute__((section(".bss.emi")));

static void root_at(unsigned phase, complex_float *root) {
    phase &= (FFT_N - 1);
    if (phase < 64) {
        root->re = k_base_cos[phase];
        root->im = k_base_sin[phase];
        return;
    }

    phase -= 64;
    root->re = -k_base_cos[phase];
    root->im = -k_base_sin[phase];
}

static void clear_matrix(complex_float *matrix) {
    unsigned i;
    for (i = 0; i < FFT_WORDS; ++i) {
        matrix[i].re = 0.0f;
        matrix[i].im = 0.0f;
    }
}

static void copy_matrix(complex_float *dst, const complex_float *src) {
    unsigned i;
    for (i = 0; i < FFT_WORDS; ++i) {
        dst[i] = src[i];
    }
}

static void set_spike(complex_float *matrix, unsigned fy, unsigned fx, float re, float im) {
    clear_matrix(matrix);
    matrix[fy * FFT_N + fx].re = re;
    matrix[fy * FFT_N + fx].im = im;
}

static unsigned next_state(unsigned state) {
    return state * 1664525u + 1013904223u;
}

static int rand_small(unsigned *state) {
    *state = next_state(*state);
    return (int)((*state >> 16) % 9u) - 4;
}

static void fill_case(int case_id, complex_float *matrix) {
    unsigned x;
    unsigned y;
    unsigned state = 1u;

    clear_matrix(matrix);

    switch (case_id) {
        case 0:
            return;
        case 1:
            matrix[0].re = 1.0f;
            return;
        case 2:
            for (y = 0; y < FFT_N; ++y) {
                for (x = 0; x < FFT_N; ++x) {
                    matrix[y * FFT_N + x].re = 1.0f;
                    matrix[y * FFT_N + x].im = 1.0f;
                }
            }
            return;
        case 3:
            for (y = 0; y < FFT_N; ++y) {
                for (x = 0; x < FFT_N; ++x) {
                    matrix[y * FFT_N + x].re = ((x + y) & 1u) ? -1.0f : 1.0f;
                }
            }
            return;
        case 4:
            for (y = 0; y < FFT_N; ++y) {
                for (x = 0; x < FFT_N; ++x) {
                    complex_float root;
                    root_at((3u * x + 5u * y) & (FFT_N - 1), &root);
                    matrix[y * FFT_N + x] = root;
                }
            }
            return;
        case 5:
            for (y = 0; y < FFT_N; ++y) {
                for (x = 0; x < FFT_N; ++x) {
                    matrix[y * FFT_N + x].re = (float)rand_small(&state);
                    matrix[y * FFT_N + x].im = (float)rand_small(&state);
                }
            }
            return;
        default:
            return;
    }
}

static void analytic_expected(int case_id, complex_float *matrix) {
    clear_matrix(matrix);

    switch (case_id) {
        case 0:
            return;
        case 1: {
            unsigned i;
            for (i = 0; i < FFT_WORDS; ++i) {
                matrix[i].re = 1.0f;
            }
            return;
        }
        case 2:
            set_spike(matrix, 0u, 0u, (float)FFT_DC_SCALE, (float)FFT_DC_SCALE);
            return;
        case 3:
            set_spike(matrix, 64u, 64u, (float)FFT_DC_SCALE, 0.0f);
            return;
        case 4:
            set_spike(matrix, 5u, 3u, (float)FFT_DC_SCALE, 0.0f);
            return;
        default:
            return;
    }
}

static void reference_dft1d(const complex_float *in, complex_float *out) {
    unsigned k;
    unsigned n;

    for (k = 0; k < FFT_N; ++k) {
        double sum_re = 0.0;
        double sum_im = 0.0;

        for (n = 0; n < FFT_N; ++n) {
            complex_float root;
            double a_re = in[n].re;
            double a_im = in[n].im;

            root_at((k * n) & (FFT_N - 1), &root);
            sum_re += a_re * root.re + a_im * root.im;
            sum_im += a_im * root.re - a_re * root.im;
        }

        out[k].re = (float)sum_re;
        out[k].im = (float)sum_im;
    }
}

static void reference_fft2d(const complex_float *in, complex_float *out) {
    unsigned row;
    unsigned col;
    unsigned ky;

    for (row = 0; row < FFT_N; ++row) {
        reference_dft1d(&in[row * FFT_N], &g_ref_scratch[row * FFT_N]);
    }

    for (col = 0; col < FFT_N; ++col) {
        complex_float column_in[FFT_N];

        for (row = 0; row < FFT_N; ++row) {
            column_in[row] = g_ref_scratch[row * FFT_N + col];
        }

        reference_dft1d(column_in, g_ref_column);

        for (ky = 0; ky < FFT_N; ++ky) {
            out[ky * FFT_N + col] = g_ref_column[ky];
        }
    }
}

static int round_nearest(float value) {
    if (value >= 0.0f) {
        return (int)(value + 0.5f);
    }
    return (int)(value - 0.5f);
}

static int max_abs_int(int value) {
    return value < 0 ? -value : value;
}

typedef struct {
    int fd;
    unsigned word;
    int bytes_left;
    int eof;
    int pushback;
    int has_pushback;
} packed_text_reader;

static int reader_fill_word(packed_text_reader *reader) {
    unsigned word = 0;
    int rc;

    if (reader->eof) {
        return 0;
    }

    rc = read(reader->fd, &word, 1);
    if (rc < 0) {
        return -1;
    }
    if (rc == 0 && word == 0) {
        reader->eof = 1;
        return 0;
    }

    reader->word = word;
    reader->bytes_left = 4;
    if (rc == 0) {
        reader->eof = 1;
    }
    return 1;
}

static int reader_getc(packed_text_reader *reader) {
    int ch;

    if (reader->has_pushback) {
        reader->has_pushback = 0;
        return reader->pushback;
    }

    while (reader->bytes_left == 0) {
        int status = reader_fill_word(reader);
        if (status <= 0) {
            return -1;
        }
    }

    ch = (int)(reader->word & 0xffu);
    reader->word >>= 8;
    --reader->bytes_left;
    if (ch == 0) {
        reader->bytes_left = 0;
        return -1;
    }
    return ch;
}

static void reader_ungetc(packed_text_reader *reader, int ch) {
    reader->has_pushback = 1;
    reader->pushback = ch;
}

static int parse_next_int_reader(packed_text_reader *reader, int *value) {
    int ch;
    int sign = 1;
    int parsed = 0;

    do {
        ch = reader_getc(reader);
    } while (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t');

    if (ch < 0) {
        return 0;
    }
    if (ch == '-') {
        sign = -1;
        ch = reader_getc(reader);
    } else if (ch == '+') {
        ch = reader_getc(reader);
    }
    if (ch < '0' || ch > '9') {
        return -1;
    }

    while (ch >= '0' && ch <= '9') {
        parsed = parsed * 10 + (ch - '0');
        ch = reader_getc(reader);
    }

    if (ch >= 0) {
        reader_ungetc(reader, ch);
    }

    *value = sign * parsed;
    return 1;
}

static int load_matrix_from_file(const char *path, complex_float *matrix) {
    packed_text_reader reader;
    unsigned i;
    int status = 0;

    reader.fd = open(path, O_RDONLY);
    if (reader.fd < 0) {
        return 1;
    }
    reader.word = 0;
    reader.bytes_left = 0;
    reader.eof = 0;
    reader.pushback = 0;
    reader.has_pushback = 0;

    for (i = 0; i < FFT_WORDS; ++i) {
        int re;
        int im;

        if (parse_next_int_reader(&reader, &re) != 1 || parse_next_int_reader(&reader, &im) != 1) {
            status = 1;
            break;
        }
        matrix[i].re = (float)re;
        matrix[i].im = (float)im;
    }

    close(reader.fd);
    return status;
}

static unsigned fnv1a_step(unsigned hash, unsigned value) {
    return (hash ^ value) * 16777619u;
}

static const char *open_test_name(int test_id) {
    switch (test_id) {
        case 1:
            return "open_tests/test_01_open";
        case 2:
            return "open_tests/test_02_open";
        default:
            return "open_tests/unknown";
    }
}

static int expected_open_hashes(int test_id, unsigned *hash_a, unsigned *hash_b) {
    switch (test_id) {
        case 1:
            *hash_a = 0xeb270348u;
            *hash_b = 0x53800634u;
            return 0;
        case 2:
            *hash_a = 0xecca1e32u;
            *hash_b = 0xc2f82c8eu;
            return 0;
        default:
            return 1;
    }
}

static int compare_expected_digest(const complex_float *actual, int test_id, const char *label) {
    unsigned i;
    unsigned expected_a;
    unsigned expected_b;
    unsigned got_a = 2166136261u;
    unsigned got_b = 33554467u;

    if (expected_open_hashes(test_id, &expected_a, &expected_b)) {
        printf("Unknown open test %d\n", test_id);
        return 1;
    }

    for (i = 0; i < FFT_WORDS; ++i) {
        unsigned idx = i;
        unsigned got_re = (unsigned)round_nearest(actual[i].re);
        unsigned got_im = (unsigned)round_nearest(actual[i].im);

        got_a = fnv1a_step(got_a, idx);
        got_a = fnv1a_step(got_a, got_re);
        got_a = fnv1a_step(got_a, got_im);

        got_b = fnv1a_step(got_b, idx ^ 0x9e3779b9u);
        got_b = fnv1a_step(got_b, got_im ^ 0x85ebca6bu);
        got_b = fnv1a_step(got_b, got_re ^ 0xc2b2ae35u);
    }

    if (got_a != expected_a || got_b != expected_b) {
        printf("%s FAIL digest_a=%08x expected_a=%08x digest_b=%08x expected_b=%08x\n",
               label,
               got_a,
               expected_a,
               got_b,
               expected_b);
        return 1;
    }

    return 0;
}

static int check_case(int case_id) {
    unsigned i;
    int changed_values = 0;
    int max_abs_diff = 0;

    fill_case(case_id, g_template);
    copy_matrix(g_matrix, g_template);

    if (case_id == 5) {
        reference_fft2d(g_template, g_expected);
    } else {
        analytic_expected(case_id, g_expected);
    }

    contest_fft_2d(g_matrix, k_base_cos, k_base_sin);

    for (i = 0; i < FFT_WORDS; ++i) {
        int got_re = round_nearest(g_matrix[i].re);
        int got_im = round_nearest(g_matrix[i].im);
        int exp_re = round_nearest(g_expected[i].re);
        int exp_im = round_nearest(g_expected[i].im);
        int diff_re = got_re - exp_re;
        int diff_im = got_im - exp_im;
        int abs_re = max_abs_int(diff_re);
        int abs_im = max_abs_int(diff_im);

        if (abs_re || abs_im) {
            ++changed_values;
            if (abs_re > max_abs_diff) {
                max_abs_diff = abs_re;
            }
            if (abs_im > max_abs_diff) {
                max_abs_diff = abs_im;
            }
        }
    }

    if (changed_values) {
        printf("FFT max_abs_diff=%d, changed_values=%d\n", max_abs_diff, changed_values);
        return 1;
    }

    return 0;
}

static int check_open_file(int test_id, const char *input_path) {
    if (load_matrix_from_file(input_path, g_matrix)) {
        printf("%s FAIL malformed input data\n", open_test_name(test_id));
        return 1;
    }

    contest_fft_2d(g_matrix, k_base_cos, k_base_sin);
    return compare_expected_digest(g_matrix, test_id, open_test_name(test_id));
}

static void usage(const char *prog) {
    printf("Usage: %s open-file <test_id> <input_path>\n", prog);
    printf("       %s check <case_id>\n", prog);
    printf("       %s bench <case_id> <iters>\n", prog);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "open-file") == 0) {
        if (argc < 4) {
            usage(argv[0]);
            return 1;
        }
        return check_open_file(atoi(argv[2]), argv[3]);
    }

    if (strcmp(argv[1], "check") == 0) {
        if (argc < 3) {
            usage(argv[0]);
            return 1;
        }
        int case_id = atoi(argv[2]);
        return check_case(case_id);
    }

    if (strcmp(argv[1], "bench") == 0) {
        int case_id;
        int iters;
        int iter;

        if (argc < 4) {
            usage(argv[0]);
            return 1;
        }

        case_id = atoi(argv[2]);
        iters = atoi(argv[3]);
        if (iters < 1) {
            return 1;
        }

        fill_case(case_id, g_template);
        for (iter = 0; iter < iters; ++iter) {
            copy_matrix(g_matrix, g_template);
            contest_fft_2d(g_matrix, k_base_cos, k_base_sin);
        }

        return 0;
    }

    usage(argv[0]);
    return 1;
}
