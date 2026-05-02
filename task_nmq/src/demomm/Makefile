NMC_PREFIX ?= /path/to/NMC-SDK-1.4.249.d4227c9.x86_64/usr/local
CC    = $(NMC_PREFIX)/bin/nmc-gcc
QEMU  = $(NMC_PREFIX)/bin/nmc-qemu

CFLAGS    = -O2 -Wall -mnmc4 -D__NM6408__
INCLUDES  = -I$(NMC_PREFIX)/nmc/include
LIBS      = -L$(NMC_PREFIX)/nmc/lib/nmc4 -lsemihost
LDS       = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../../ENV/nmq_qemu_linux_bundle_participant/example/qemu-nm6408.lds
LINKFLAGS = -Wl,-T$(LDS)

TARGET  = solution.abs
SOURCES = harness_argv.c solution.S

BENCH_ITERS ?= 1000
BENCH_M ?= 24
BENCH_K ?= 24
BENCH_N ?= 24

# sample-1: 2×3 × 3×2  →  C = [[58,64],[139,154]]
S1_ARGS = 2 3 2  1 2 3 4 5 6  7 8 9 10 11 12

# sample-2: 3×3 × 3×3  →  C = [[-2,5,4],[19,3,17],[20,-7,10]]
S2_ARGS = 3 3 3  2 0 -1 3 1 4 0 -2 5  1 2 3 0 1 0 4 -1 2

.PHONY: all clean run test test1 test2 verify-large bench dump

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile
	$(CC) -o $@ $(SOURCES) $(CFLAGS) $(INCLUDES) $(LIBS) $(LINKFLAGS)
	@echo "Build successful: $@"

test1: $(TARGET)
	@printf '%s\n' '58 64' '139 154' > /tmp/nmc_expected.txt; \
	$(QEMU) $(TARGET) $(S1_ARGS) 2>/dev/null | grep -v 'Successful end' | grep . > /tmp/nmc_got.txt; \
	if diff -q /tmp/nmc_expected.txt /tmp/nmc_got.txt >/dev/null 2>&1; then \
		echo "PASS sample-1"; \
	else \
		echo "FAIL sample-1"; \
		echo "  expected:"; cat /tmp/nmc_expected.txt; \
		echo "  got:";      cat /tmp/nmc_got.txt; exit 1; \
	fi

test2: $(TARGET)
	@printf '%s\n' '-2 5 4' '19 3 17' '20 -7 10' > /tmp/nmc_expected.txt; \
	$(QEMU) $(TARGET) $(S2_ARGS) 2>/dev/null | grep -v 'Successful end' | grep . > /tmp/nmc_got.txt; \
	if diff -q /tmp/nmc_expected.txt /tmp/nmc_got.txt >/dev/null 2>&1; then \
		echo "PASS sample-2"; \
	else \
		echo "FAIL sample-2"; \
		echo "  expected:"; cat /tmp/nmc_expected.txt; \
		echo "  got:";      cat /tmp/nmc_got.txt; exit 1; \
	fi

test: test1 test2

verify-large: $(TARGET)
	@for dims in \
		"24 24 24" "32 32 32" "48 48 48" "64 64 64" "128 128 128" \
		"24 32 40" "32 24 48" "48 64 32" "64 48 96" "64 64 128" "128 64 64"; do \
		set -- $$dims; \
		if $(QEMU) $(TARGET) checkgen $$1 $$2 $$3 >/dev/null 2>&1; then \
			echo "PASS $$1x$$2x$$3"; \
		else \
			echo "FAIL $$1x$$2x$$3"; exit 1; \
		fi; \
	done

run: $(TARGET)
	$(QEMU) $(TARGET) $(S1_ARGS)

# Bench: host wall-clock time across multiple matrix shapes.
# Override BENCH_ITERS to change iteration count per shape.
# Single QEMU invocation covers all shapes (multibench mode) — eliminates
# per-shape startup overhead while preserving the same workload.
BENCH_SHAPES = 2 2 2  4 4 4  8 8 8  12 12 12  16 16 16  24 24 24 \
               2 24 2  4 16 4  2 32 2  6 16 6

bench: $(TARGET)
	@echo "Benchmarking $(BENCH_ITERS) iters/shape (single QEMU invocation)..."
	@printf '%-14s %10s\n' "shapes" "total_ms"
	@printf '%-14s %10s\n' "------" "--------"
	@START=$$(date +%s%N); \
	$(QEMU) $(TARGET) multi $(BENCH_ITERS) $(BENCH_SHAPES) >/dev/null 2>&1; \
	END=$$(date +%s%N); \
	TOTAL_NS=$$((END - START)); \
	printf '%-14s %10s\n' "all" "$$(( TOTAL_NS / 1000000 ))"

clean:
	rm -f *.abs *.map *.dump

dump: $(TARGET)
	$(NMC_PREFIX)/bin/nmc-objdump -d -x -s -l $(TARGET) > dump.dump
