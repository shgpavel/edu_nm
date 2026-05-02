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

BENCH_ITERS ?= 2500
BENCH_BASE ?= 4294967291
BENCH_EXP ?= 4294967295
BENCH_MOD ?= 4294967279

.PHONY: all clean run test bench dump

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile
	$(CC) -o $@ $(SOURCES) $(CFLAGS) $(INCLUDES) $(LIBS) $(LINKFLAGS)
	@echo "Build successful: $@"

test: $(TARGET)
	@PASS=0; FAIL=0; \
	for case_data in \
		"2 10 1000 24" \
		"0 0 7 1" \
		"123456789 0 1 0" \
		"4294967295 1 4294967294 1" \
		"2 32 1000000007 294967268" \
		"3 0 5 1" \
		"7 1 13 7" \
		"2 0 1000000007 1" \
	; do \
		set -- $$case_data; \
		BASE=$$1; EXP=$$2; MOD=$$3; EXPECTED=$$4; \
		RESULT=$$($(QEMU) $(TARGET) $$BASE $$EXP $$MOD 2>/dev/null | head -1); \
		if [ "$$RESULT" = "$$EXPECTED" ]; then \
			echo "PASS: pow_mod($$BASE, $$EXP, $$MOD) = $$RESULT"; PASS=$$((PASS+1)); \
		else \
			echo "FAIL: pow_mod($$BASE, $$EXP, $$MOD), expected $$EXPECTED, got $$RESULT"; FAIL=$$((FAIL+1)); \
		fi; \
	done; \
	echo "Results: $$PASS passed, $$FAIL failed"; \
	[ $$FAIL -eq 0 ]

run: $(TARGET)
	$(QEMU) $(TARGET) 2 10 1000

bench: $(TARGET)
	@echo "Benchmarking $(BENCH_ITERS) iterations (base=$(BENCH_BASE), exp=$(BENCH_EXP), mod=$(BENCH_MOD))..."
	@START=$$(date +%s%N); \
	$(QEMU) $(TARGET) $(BENCH_BASE) $(BENCH_EXP) $(BENCH_MOD) $(BENCH_ITERS) >/dev/null 2>&1; \
	END=$$(date +%s%N); \
	TOTAL_NS=$$((END - START)); \
	NS_PER=$$(( TOTAL_NS / $(BENCH_ITERS) )); \
	echo "total wall time: $$(( TOTAL_NS / 1000000 )) ms"; \
	echo "ns/call (QEMU proxy): $$NS_PER ns"

clean:
	rm -f *.abs *.map *.dump

dump: $(TARGET)
	$(NMC_PREFIX)/bin/nmc-objdump -d -x -s -l $(TARGET) > dump.dump
