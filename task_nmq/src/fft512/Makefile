NMC_PREFIX ?= /path/to/NMC-SDK-1.4.249.d4227c9.x86_64/usr/local
CC    = $(NMC_PREFIX)/bin/nmc-gcc
QEMU  = $(NMC_PREFIX)/bin/nmc-qemu

CFLAGS    = -O2 -Wall -mnmc4 -D__NM6408__
INCLUDES  = -I$(NMC_PREFIX)/nmc/include
LIBS      = -L$(NMC_PREFIX)/nmc/lib/nmc4 -lsemihost
LDS       = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../../ENV/nmq_qemu_linux_bundle_participant/example/qemu-nm6408.lds
LINKFLAGS = -Wl,-T$(LDS)

TARGET  = solution.abs
TARGET_SOURCES = harness_main.c solution.S
OPEN_TESTS = test_01_open test_02_open

BENCH_ITERS ?= 1
BENCH_CASE  ?= 5

.PHONY: all clean run test bench dump

all: $(TARGET)

# Hand-maintained during optimization; keep solution.S stable.

$(TARGET): $(TARGET_SOURCES) Makefile
	$(CC) -o $@ $(TARGET_SOURCES) $(CFLAGS) $(INCLUDES) $(LIBS) $(LINKFLAGS)
	@echo "Build successful: $@"

test: $(TARGET)
	@for test_name in $(OPEN_TESTS); do \
		test_dir="open_tests/$$test_name"; \
		if $(QEMU) $(TARGET) open-file "$$test_dir/input.txt" "$$test_dir/expected_output.txt" >/tmp/fft512_test.txt 2>&1; then \
			echo "PASS $$test_dir"; \
		else \
			cat /tmp/fft512_test.txt; \
			exit 1; \
		fi; \
	done

run: $(TARGET)
	$(QEMU) $(TARGET) check 4

bench: $(TARGET)
	@echo "Benchmarking $(BENCH_ITERS) iterations (case=$(BENCH_CASE))..."
	@START=$$(date +%s%N); \
	$(QEMU) $(TARGET) bench $(BENCH_CASE) $(BENCH_ITERS) >/dev/null 2>&1; \
	END=$$(date +%s%N); \
	TOTAL_NS=$$((END - START)); \
	NS_PER=$$(( TOTAL_NS / $(BENCH_ITERS) )); \
	echo "total wall time: $$(( TOTAL_NS / 1000000 )) ms"; \
	echo "ns/call (QEMU proxy): $$NS_PER ns"

clean:
	rm -f *.abs *.map *.dump *.o

dump: $(TARGET)
	$(NMC_PREFIX)/bin/nmc-objdump -d -x -s -l $(TARGET) > dump.dump
