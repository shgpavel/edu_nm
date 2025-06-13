CFLAGS  := -Wall -Wextra -O2 -std=c23 -march=native $(CFLAGS)
LDFLAGS := -lc -lm -flto -fuse-ld=lld $(LDFLAGS)

LDFLAGS += -lmkl_core -ldl -lmkl_intel_lp64 -ltbb \
	-lmkl_tbb_thread -lpthread

CC = clang
LD = clang

LD_LIBRARY_PATH := -L/usr/lib
INCL_PATH       := -Iinclude/

OUT_O_DIR := build
SRCS_DIRS := . methods types
SRCS      := $(wildcard $(addsuffix /*.c, $(SRCS_DIRS)))
OBJS      := $(addprefix $(OUT_O_DIR)/,$(SRCS:.c=.o))
DEPS      := $(OBJS:.o=.d)

ifdef DEBUG
	CFLAGS += -DDEBUG -O1 -g
endif

ifdef SANITIZE
	CFLAGS  += -fsanitize=address -fsanitize=undefined
	LDFLAGS += -fsanitize=address -fsanitize=undefined
endif

ifdef PGO_GEN
	CFLAGS += -fprofile-generate
	LDFLAGS += -fprofile-generate
endif

ifdef PGO_USE
	CFLAGS += -fprofile-use=$(PGO_DATA)
	LDFLAGS += -fprofile-use=$(PGO_DATA)
endif

.PHONY: all clean-build clean-all clean \
	pgo pgo-gen pgo-run pgo-merge pgo-use

all: $(OUT_O_DIR)/main

$(OUT_O_DIR)/main: $(OBJS)
	$(LD) $^ -o $@ $(LDFLAGS) $(LD_LIBRARY_PATH)

$(OBJS): $(OUT_O_DIR)/%.o: %.c Makefile
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(INCL_PATH) -c $< -o $@

$(DEPS): $(OUT_O_DIR)/%.d: %.c Makefile
	@mkdir -p $(@D)
	$(CC) -E $(CFLAGS) $(INCL_PATH) $< -MM -MT $(@:.d=.o) > $@

clean: clean-all

clean-all: clean-build
	rm -f $(PGO_RAW) $(PGO_DATA)

clean-build:
	rm -rf $(OBJS) $(DEPS) $(OUT_O_DIR)

NODEPS = clean clean-all clean-build

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
include $(DEPS)
endif

PGO_RAW  = default.profraw
PGO_DATA = default.profdata

pgo: clean pgo-gen pgo-run pgo-merge pgo-use
pgo-gen:
	$(MAKE) all PGO_GEN=1

pgo-run:
	@env > /tmp/original_env && \
	(env -i bash -c \
	". /opt/intel/oneapi/setvars.sh && LLVM_PROFILE_FILE=$(PGO_RAW) ./build/main") \
	&& eval "$$(cat /tmp/original_env)"

pgo-merge:
	llvm-profdata merge -sparse $(PGO_RAW) -o $(PGO_DATA)

pgo-use:
	$(MAKE) clean-build
	$(MAKE) all PGO_USE=1 PGO_DATA=$(PGO_DATA)
