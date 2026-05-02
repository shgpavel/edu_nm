NMC_PREFIX ?= /path/to/NMC-SDK-1.4.249.d4227c9.x86_64/usr/local
CC    = $(NMC_PREFIX)/bin/nmc-gcc
QEMU  = $(NMC_PREFIX)/bin/nmc-qemu

CFLAGS    = -O2 -Wall -mnmc4 -D__NM6408__
INCLUDES  = -I$(NMC_PREFIX)/nmc/include
LIBS      = -L$(NMC_PREFIX)/nmc/lib/nmc4 -lsemihost
LDS       = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../../ENV/nmq_qemu_linux_bundle_participant/example/qemu-nm6408.lds
LINKFLAGS = -Wl,-T$(LDS)

TARGET  = solution.abs
SOURCES = harness_main.S solution.S

BENCH_ITERS ?= 1000
BENCH_N ?= 12

.PHONY: all clean run test bench dump

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile
	$(CC) -o $@ $(SOURCES) $(CFLAGS) $(INCLUDES) $(LIBS) $(LINKFLAGS)
	@echo "Build successful: $@"

test: $(TARGET)
	@CASE1_EXPECTED="7.111111 -3.222222"; \
	CASE1_RESULT=$$($(QEMU) $(TARGET) 2 2 1 5 7 11 13 2>/dev/null | head -1); \
	if [ "$$CASE1_RESULT" = "$$CASE1_EXPECTED" ]; then \
		echo "PASS case-1"; \
	else \
		echo "FAIL case-1: expected $$CASE1_EXPECTED, got $$CASE1_RESULT"; \
		exit 1; \
	fi
	@CASE2_EXPECTED="1.000000 2.000000 3.000000"; \
	CASE2_RESULT=$$($(QEMU) $(TARGET) 3 1 1 1 2 5 5 2 3 8 6 27 32 2>/dev/null | head -1); \
	if [ "$$CASE2_RESULT" = "$$CASE2_EXPECTED" ]; then \
		echo "PASS case-2"; \
	else \
		echo "FAIL case-2: expected $$CASE2_EXPECTED, got $$CASE2_RESULT"; \
		exit 1; \
	fi
	@CASE3_EXPECTED="SINGULAR"; \
	CASE3_RESULT=$$($(QEMU) $(TARGET) 2 1 2 2 4 3 6 2>/dev/null | head -1); \
	if [ "$$CASE3_RESULT" = "$$CASE3_EXPECTED" ]; then \
		echo "PASS case-3"; \
	else \
		echo "FAIL case-3: expected $$CASE3_EXPECTED, got $$CASE3_RESULT"; \
		exit 1; \
	fi

run: $(TARGET)
	$(QEMU) $(TARGET) 2 2 1 5 7 11 13

bench: $(TARGET)
	@echo "Benchmarking $(BENCH_ITERS) iterations (generated n=$(BENCH_N))..."
	@START=$$(date +%s%N); \
	$(QEMU) $(TARGET) gen $(BENCH_N) $(BENCH_ITERS) >/dev/null 2>&1; \
	END=$$(date +%s%N); \
	TOTAL_NS=$$((END - START)); \
	NS_PER=$$(( TOTAL_NS / $(BENCH_ITERS) )); \
	echo "total wall time: $$(( TOTAL_NS / 1000000 )) ms"; \
	echo "ns/call (QEMU proxy): $$NS_PER ns"

clean:
	rm -f *.abs *.map *.dump

dump: $(TARGET)
	$(NMC_PREFIX)/bin/nmc-objdump -d -x -s -l $(TARGET) > dump.dump
