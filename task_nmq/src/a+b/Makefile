NMC_PREFIX ?= /path/to/NMC-SDK-1.4.249.d4227c9.x86_64/usr/local
CC    = $(NMC_PREFIX)/bin/nmc-gcc
QEMU  = $(NMC_PREFIX)/bin/nmc-qemu

CFLAGS   = -O2 -Wall -mnmc4 -D__NM6408__
INCLUDES = -I$(NMC_PREFIX)/nmc/include
LIBS     = -L$(NMC_PREFIX)/nmc/lib/nmc4 -lsemihost
LDS      = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../../ENV/nmq_qemu_linux_bundle_participant/example/qemu-nm6408.lds
LINKFLAGS = -Wl,-T$(LDS)

TARGET  = solution.abs
SOURCES = harness_main.S solution.S

BENCH_ITERS ?= 3000
BENCH_A ?= 2147483640
BENCH_B ?= -2147480000

.PHONY: all clean run test bench

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile
	$(CC) -o $@ $(SOURCES) $(CFLAGS) $(INCLUDES) $(LIBS) $(LINKFLAGS)
	@echo "Build successful: $@"

test: $(TARGET)
	@PASS=0; FAIL=0; \
	for case_data in \
		"1 2 3" \
		"0 0 0" \
		"-1 1 0" \
		"-5 -3 -8" \
		"2147483647 0 2147483647" \
		"-2147483648 0 -2147483648" \
		"2147483640 -2147480000 3640" \
		"1000000 -1000000 0" \
	; do \
		set -- $$case_data; A=$$1; B=$$2; EXPECTED=$$3; \
		RESULT=$$($(QEMU) $(TARGET) $$A $$B 2>/dev/null | head -1); \
		if [ "$$RESULT" = "$$EXPECTED" ]; then \
			echo "PASS: a_plus_b($$A, $$B) = $$RESULT"; PASS=$$((PASS+1)); \
		else \
			echo "FAIL: a_plus_b($$A, $$B), expected $$EXPECTED, got $$RESULT"; FAIL=$$((FAIL+1)); \
		fi; \
	done; \
	echo "Results: $$PASS passed, $$FAIL failed"; \
	[ $$FAIL -eq 0 ]

run: $(TARGET)
	$(QEMU) $(TARGET) 1 2

# bench: time BENCH_ITERS calls from the host (clock() is a stub inside QEMU).
# Reports wall-time ns/call as a relative performance proxy.
bench: $(TARGET)
	@echo "Benchmarking $(BENCH_ITERS) iterations (a=$(BENCH_A), b=$(BENCH_B))..."
	@START=$$(date +%s%N); \
	$(QEMU) $(TARGET) $(BENCH_A) $(BENCH_B) $(BENCH_ITERS) >/dev/null 2>&1; \
	END=$$(date +%s%N); \
	TOTAL_NS=$$((END - START)); \
	NS_PER=$$(( TOTAL_NS / $(BENCH_ITERS) )); \
	echo "total wall time: $$(( TOTAL_NS / 1000000 )) ms"; \
	echo "ns/call (QEMU proxy): $$NS_PER ns"

clean:
	rm -f *.abs *.map *.dump

dump: $(TARGET)
	$(NMC_PREFIX)/bin/nmc-objdump -d -x -s -l $(TARGET) > dump.dump
