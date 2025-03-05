CFLAGS := -Wall -Wextra -O3 -std=c23 -march=native -fsanitize=address -fsanitize=undefined $(CFLAGS)
LDFLAGS := -lc -lm -fuse-ld=mold -fsanitize=address -fsanitize=undefined $(LDFLAGS)

# 

CC = clang
LD = clang

LD_LIBRARY_PATH = -L/usr/lib
INCL_PATH = -Iinclude/

ifeq ($(DEBUG_INFO),1)
	CFLAGS += -DDEBUG_INFO
endif

OUT_O_DIR = build
SRCS_DIRS = . methods
SRCS = $(wildcard $(addsuffix /*.c, $(SRCS_DIRS)))
OBJS := $(addprefix $(OUT_O_DIR)/,$(SRCS:.c=.o))
DEPS = $(OBJS:.o=.d)

.PHONY: all
all: $(OUT_O_DIR)/main

$(OUT_O_DIR)/main: $(OBJS)
	$(LD) $^ -o $@ $(LDFLAGS) $(LD_LIBRARY_PATH)

$(OBJS): $(OUT_O_DIR)/%.o: %.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(INCL_PATH) -c $< -o $@

$(DEPS): $(OUT_O_DIR)/%.d: %.c
	@mkdir -p $(@D)
	$(CC) -E $(CFLAGS) $(INCL_PATH) $< -MM -MT $(@:.d=.o) > $@

.PHONY: clean
clean:
	rm -rf $(OBJS) $(DEPS) $(OUT_O_DIR)

NODEPS = clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
include $(DEPS)
endif
