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

BENCH_ITERS ?= 3000
BENCH_N ?= 100000

.PHONY: all clean run test bench dump

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile
	$(CC) -o $@ $(SOURCES) $(CFLAGS) $(INCLUDES) $(LIBS) $(LINKFLAGS)
	@echo "Build successful: $@"

test: $(TARGET)
	@PASS=0; FAIL=0; \
	for case_data in \
		"0 1" \
		"1 1" \
		"2 2" \
		"5 2" \
		"10 8" \
		"20 4" \
		"25 4" \
		"100 4" \
		"1000 2" \
		"10000 8" \
	; do \
		set -- $$case_data; N=$$1; EXPECTED=$$2; \
		RESULT=$$($(QEMU) $(TARGET) $$N 2>/dev/null | head -1); \
		if [ "$$RESULT" = "$$EXPECTED" ]; then \
			echo "PASS: lastfct($$N) = $$RESULT"; PASS=$$((PASS+1)); \
		else \
			echo "FAIL: lastfct($$N), expected $$EXPECTED, got $$RESULT"; FAIL=$$((FAIL+1)); \
		fi; \
	done; \
	echo "Results: $$PASS passed, $$FAIL failed"; \
	[ $$FAIL -eq 0 ]

run: $(TARGET)
	$(QEMU) $(TARGET) 20

bench: $(TARGET)
	@echo "Benchmarking $(BENCH_ITERS) iterations (n=$(BENCH_N))..."
	@START=$$(date +%s%N); \
	$(QEMU) $(TARGET) $(BENCH_N) $(BENCH_ITERS) >/dev/null 2>&1; \
	END=$$(date +%s%N); \
	TOTAL_NS=$$((END - START)); \
	NS_PER=$$(( TOTAL_NS / $(BENCH_ITERS) )); \
	echo "total wall time: $$(( TOTAL_NS / 1000000 )) ms"; \
	echo "ns/call (QEMU proxy): $$NS_PER ns"

clean:
	rm -f *.abs *.map *.dump

dump: $(TARGET)
	$(NMC_PREFIX)/bin/nmc-objdump -d -x -s -l $(TARGET) > dump.dump
