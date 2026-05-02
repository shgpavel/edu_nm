NMC_PREFIX ?= /path/to/NMC-SDK-1.4.249.d4227c9.x86_64/usr/local
CC    = $(NMC_PREFIX)/bin/nmc-gcc
QEMU  = $(NMC_PREFIX)/bin/nmc-qemu

CFLAGS   = -O2 -Wall -mnmc4 -D__NM6408__
INCLUDES = -I$(NMC_PREFIX)/nmc/include
LIBS     = -L$(NMC_PREFIX)/nmc/lib/nmc4 -lsemihost
LDS      = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../../ENV/nmq_qemu_linux_bundle_participant/example/qemu-nm6408.lds
LINKFLAGS = -Wl,-T$(LDS)

TARGET  = solution.abs
SOURCES = harness_argv.S solution.S

BENCH_ITERS ?= 1000

.PHONY: all clean run test bench

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile
	$(CC) -o $@ $(SOURCES) $(CFLAGS) $(INCLUDES) $(LIBS) $(LINKFLAGS)
	@echo "Build successful: $@"

# RC4 test vectors — format: key_len k0..kN data_len d0..dM -> expected output bytes
# Vector 1: RFC 6229 key=1..5 (0x01..0x05), plaintext=0 0 0 0 0
# Vector 2: key=all zeros (1 byte), data=0 0
# Vector 3: key="Key" (75 101 121), data="Plaintext" (80 108 97 105 110 116 101 120 116)
test: $(TARGET)
	@PASS=0; FAIL=0; \
	V1_RESULT=$$($(QEMU) $(TARGET) 5 1 2 3 4 5 5 0 0 0 0 0 2>/dev/null | head -1 | tr -s ' ' | sed 's/ *$$//'); \
	V1_EXPECTED="178 57 99 5 240"; \
	if [ "$$V1_RESULT" = "$$V1_EXPECTED" ]; then \
		echo "PASS: rc4 vector-1"; PASS=$$((PASS+1)); \
	else \
		echo "FAIL: rc4 vector-1, expected [$$V1_EXPECTED], got [$$V1_RESULT]"; FAIL=$$((FAIL+1)); \
	fi; \
	V2_RESULT=$$($(QEMU) $(TARGET) 1 0 2 0 0 2>/dev/null | head -1 | tr -s ' ' | sed 's/ *$$//'); \
	V2_EXPECTED="222 24"; \
	if [ "$$V2_RESULT" = "$$V2_EXPECTED" ]; then \
		echo "PASS: rc4 vector-2"; PASS=$$((PASS+1)); \
	else \
		echo "FAIL: rc4 vector-2, expected [$$V2_EXPECTED], got [$$V2_RESULT]"; FAIL=$$((FAIL+1)); \
	fi; \
	V3_RESULT=$$($(QEMU) $(TARGET) 3 75 101 121 9 80 108 97 105 110 116 101 120 116 2>/dev/null | head -1 | tr -s ' ' | sed 's/ *$$//'); \
	V3_EXPECTED="187 243 22 232 217 64 175 10 211"; \
	if [ "$$V3_RESULT" = "$$V3_EXPECTED" ]; then \
		echo "PASS: rc4 vector-3"; PASS=$$((PASS+1)); \
	else \
		echo "FAIL: rc4 vector-3, expected [$$V3_EXPECTED], got [$$V3_RESULT]"; FAIL=$$((FAIL+1)); \
	fi; \
	echo "Results: $$PASS passed, $$FAIL failed"; \
	[ $$FAIL -eq 0 ]

bench: $(TARGET)
	@echo "=== RC4 benchmark (iters=$(BENCH_ITERS)) ==="; \
	median3() { printf '%s\n' "$$1" "$$2" "$$3" | sort -n | sed -n '2p'; }; \
	run_once() { S=$$(date +%s%N); $(QEMU) $(TARGET) $$@ $(BENCH_ITERS) >/dev/null 2>/dev/null; E=$$(date +%s%N); echo $$(( (E - S) / $(BENCH_ITERS) )); }; \
	echo "--- small (key=16, data=16) ---"; \
	S1=$$(run_once gen 16 16); S2=$$(run_once gen 16 16); S3=$$(run_once gen 16 16); \
	SM=$$(median3 $$S1 $$S2 $$S3); \
	echo "  trials: $$S1 $$S2 $$S3 ns/call"; \
	echo "  median: $$SM ns/call  |  ns/byte: $$(( SM / 16 ))"; \
	echo "--- large (key=16, data=4096) ---"; \
	L1=$$(run_once gen 16 4096); L2=$$(run_once gen 16 4096); L3=$$(run_once gen 16 4096); \
	LM=$$(median3 $$L1 $$L2 $$L3); \
	echo "  trials: $$L1 $$L2 $$L3 ns/call"; \
	echo "  median: $$LM ns/call  |  ns/byte: $$(( LM / 4096 ))"; \
	echo "BENCH_NS_CALL=$$LM"

clean:
	rm -f *.abs *.map *.dump

dump: $(TARGET)
	$(NMC_PREFIX)/bin/nmc-objdump -d -x -s -l $(TARGET) > dump.dump
