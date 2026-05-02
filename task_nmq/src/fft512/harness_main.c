#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define FFT_N 512u
#define FFT_WORDS (FFT_N * FFT_N)
#define FFT_DC_SCALE (FFT_N * FFT_N)

typedef struct {
    float re;
    float im;
} complex_float;

typedef struct {
    int fd;
    unsigned word;
    int bytes_left;
    int eof;
    int pushback;
    int has_pushback;
} packed_text_reader;

void contest_fft_2d(complex_float *matrix, const float *twiddle_cos, const float *twiddle_sin);

static const double k_root_step_re = 0x1.fff62169b92dbp-1;
static const double k_root_step_im = 0x1.921d1fcdec784p-7;

static complex_float g_matrix[FFT_WORDS] __attribute__((section(".bss.emi")));
static complex_float g_template[FFT_WORDS] __attribute__((section(".bss.emi")));
static complex_float g_roots[FFT_N] __attribute__((section(".bss.emi")));
static float g_twiddle_cos[FFT_N / 2u] __attribute__((section(".bss.emi")));
static float g_twiddle_sin[FFT_N / 2u] __attribute__((section(".bss.emi")));
static int g_tables_ready __attribute__((section(".bss.emi")));

static void init_tables(void) {
    unsigned phase;
    double root_re = 1.0;
    double root_im = 0.0;

    if (g_tables_ready) {
        return;
    }

    g_roots[0].re = 1.0f;
    g_roots[0].im = 0.0f;
    g_twiddle_cos[0] = 1.0f;
    g_twiddle_sin[0] = 0.0f;

    for (phase = 1u; phase < FFT_N / 2u; ++phase) {
        double next_re = root_re * k_root_step_re - root_im * k_root_step_im;
        double next_im = root_im * k_root_step_re + root_re * k_root_step_im;

        root_re = next_re;
        root_im = next_im;
        g_roots[phase].re = (float)root_re;
        g_roots[phase].im = (float)root_im;
        g_twiddle_cos[phase] = (float)root_re;
        g_twiddle_sin[phase] = (float)root_im;
    }

    g_roots[FFT_N / 2u].re = -1.0f;
    g_roots[FFT_N / 2u].im = 0.0f;

    for (phase = 1u; phase < FFT_N / 2u; ++phase) {
        g_roots[phase + FFT_N / 2u].re = -g_roots[phase].re;
        g_roots[phase + FFT_N / 2u].im = -g_roots[phase].im;
    }

    g_tables_ready = 1;
}

static void root_at(unsigned phase, complex_float *root) {
    phase &= (FFT_N - 1u);
    *root = g_roots[phase];
}

static void clear_matrix(complex_float *matrix) {
    unsigned i;

    for (i = 0u; i < FFT_WORDS; ++i) {
        matrix[i].re = 0.0f;
        matrix[i].im = 0.0f;
    }
}

static void copy_matrix(complex_float *dst, const complex_float *src) {
    unsigned i;

    for (i = 0u; i < FFT_WORDS; ++i) {
        dst[i] = src[i];
    }
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
            for (y = 0u; y < FFT_N; ++y) {
                for (x = 0u; x < FFT_N; ++x) {
                    matrix[y * FFT_N + x].re = 1.0f;
                    matrix[y * FFT_N + x].im = 1.0f;
                }
            }
            return;
        case 3:
            for (y = 0u; y < FFT_N; ++y) {
                for (x = 0u; x < FFT_N; ++x) {
                    matrix[y * FFT_N + x].re = ((x + y) & 1u) ? -1.0f : 1.0f;
                }
            }
            return;
        case 4:
            for (y = 0u; y < FFT_N; ++y) {
                for (x = 0u; x < FFT_N; ++x) {
                    complex_float root;

                    root_at((3u * x + 5u * y) & (FFT_N - 1u), &root);
                    matrix[y * FFT_N + x] = root;
                }
            }
            return;
        case 5:
            for (y = 0u; y < FFT_N; ++y) {
                for (x = 0u; x < FFT_N; ++x) {
                    matrix[y * FFT_N + x].re = (float)rand_small(&state);
                    matrix[y * FFT_N + x].im = (float)rand_small(&state);
                }
            }
            return;
        default:
            return;
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

static void analytic_expected_at(int case_id, unsigned y, unsigned x, int *re, int *im) {
    *re = 0;
    *im = 0;

    switch (case_id) {
        case 1:
            *re = 1;
            return;
        case 2:
            if (y == 0u && x == 0u) {
                *re = (int)FFT_DC_SCALE;
                *im = (int)FFT_DC_SCALE;
            }
            return;
        case 3:
            if (y == FFT_N / 2u && x == FFT_N / 2u) {
                *re = (int)FFT_DC_SCALE;
            }
            return;
        case 4:
            if (y == 5u && x == 3u) {
                *re = (int)FFT_DC_SCALE;
            }
            return;
        default:
            return;
    }
}

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

    for (i = 0u; i < FFT_WORDS; ++i) {
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

static int compare_expected_file(const complex_float *actual, const char *expected_path, const char *label) {
    packed_text_reader reader;
    unsigned i;
    int mismatches = 0;
    int max_abs_diff = 0;
    unsigned first_idx = 0u;
    int first_got_re = 0;
    int first_got_im = 0;
    int first_exp_re = 0;
    int first_exp_im = 0;

    reader.fd = open(expected_path, O_RDONLY);
    if (reader.fd < 0) {
        printf("%s FAIL cannot open expected output\n", label);
        return 1;
    }
    reader.word = 0;
    reader.bytes_left = 0;
    reader.eof = 0;
    reader.pushback = 0;
    reader.has_pushback = 0;

    for (i = 0u; i < FFT_WORDS; ++i) {
        int exp_re;
        int exp_im;
        int got_re;
        int got_im;
        int diff_re;
        int diff_im;
        int abs_re;
        int abs_im;

        if (parse_next_int_reader(&reader, &exp_re) != 1 || parse_next_int_reader(&reader, &exp_im) != 1) {
            close(reader.fd);
            printf("%s FAIL malformed expected output\n", label);
            return 1;
        }

        got_re = round_nearest(actual[i].re);
        got_im = round_nearest(actual[i].im);
        diff_re = got_re - exp_re;
        diff_im = got_im - exp_im;
        abs_re = max_abs_int(diff_re);
        abs_im = max_abs_int(diff_im);

        if (abs_re || abs_im) {
            if (mismatches == 0) {
                first_idx = i;
                first_got_re = got_re;
                first_got_im = got_im;
                first_exp_re = exp_re;
                first_exp_im = exp_im;
            }
            ++mismatches;
            if (abs_re > max_abs_diff) {
                max_abs_diff = abs_re;
            }
            if (abs_im > max_abs_diff) {
                max_abs_diff = abs_im;
            }
        }
    }

    close(reader.fd);

    if (mismatches) {
        printf("%s FAIL first=(%u,%u) got=%d %d expected=%d %d mismatches=%d max_abs_diff=%d\n",
               label,
               first_idx / FFT_N,
               first_idx % FFT_N,
               first_got_re,
               first_got_im,
               first_exp_re,
               first_exp_im,
               mismatches,
               max_abs_diff);
        return 1;
    }

    return 0;
}

static int check_case(int case_id) {
    unsigned i;
    int mismatches = 0;
    int max_abs_diff = 0;
    unsigned first_idx = 0u;
    int first_got_re = 0;
    int first_got_im = 0;
    int first_exp_re = 0;
    int first_exp_im = 0;

    if (case_id < 0 || case_id > 4) {
        printf("Unsupported analytic case %d\n", case_id);
        return 1;
    }

    fill_case(case_id, g_template);
    copy_matrix(g_matrix, g_template);
    contest_fft_2d(g_matrix, g_twiddle_cos, g_twiddle_sin);

    for (i = 0u; i < FFT_WORDS; ++i) {
        int got_re = round_nearest(g_matrix[i].re);
        int got_im = round_nearest(g_matrix[i].im);
        int exp_re;
        int exp_im;
        int diff_re;
        int diff_im;
        int abs_re;
        int abs_im;

        analytic_expected_at(case_id, i / FFT_N, i % FFT_N, &exp_re, &exp_im);
        diff_re = got_re - exp_re;
        diff_im = got_im - exp_im;
        abs_re = max_abs_int(diff_re);
        abs_im = max_abs_int(diff_im);

        if (abs_re || abs_im) {
            if (mismatches == 0) {
                first_idx = i;
                first_got_re = got_re;
                first_got_im = got_im;
                first_exp_re = exp_re;
                first_exp_im = exp_im;
            }
            ++mismatches;
            if (abs_re > max_abs_diff) {
                max_abs_diff = abs_re;
            }
            if (abs_im > max_abs_diff) {
                max_abs_diff = abs_im;
            }
        }
    }

    if (mismatches) {
        printf("check %d FAIL first=(%u,%u) got=%d %d expected=%d %d mismatches=%d max_abs_diff=%d\n",
               case_id,
               first_idx / FFT_N,
               first_idx % FFT_N,
               first_got_re,
               first_got_im,
               first_exp_re,
               first_exp_im,
               mismatches,
               max_abs_diff);
        return 1;
    }

    return 0;
}

static int check_open_file(const char *input_path, const char *expected_path) {
    if (load_matrix_from_file(input_path, g_matrix)) {
        printf("open-file FAIL malformed input data\n");
        return 1;
    }

    contest_fft_2d(g_matrix, g_twiddle_cos, g_twiddle_sin);
    return compare_expected_file(g_matrix, expected_path, input_path);
}

static void usage(const char *prog) {
    printf("Usage: %s open-file <input_path> <expected_path>\n", prog);
    printf("       %s check <case_id>\n", prog);
    printf("       %s bench <case_id> <iters>\n", prog);
}

int main(int argc, char **argv) {
    init_tables();

    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "open-file") == 0) {
        if (argc < 4) {
            usage(argv[0]);
            return 1;
        }
        return check_open_file(argv[2], argv[3]);
    }

    if (strcmp(argv[1], "check") == 0) {
        if (argc < 3) {
            usage(argv[0]);
            return 1;
        }
        return check_case(atoi(argv[2]));
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
            contest_fft_2d(g_matrix, g_twiddle_cos, g_twiddle_sin);
        }

        return 0;
    }

    usage(argv[0]);
    return 1;
}
