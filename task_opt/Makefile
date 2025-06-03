CFLAGS  := -Wall -Wextra -O2 -std=c23 $(CFLAGS)
LDFLAGS := -lc -lm -fuse-ld=lld $(LDFLAGS)

# -fsanitize=address -fsanitize=undefined 

CC = clang
LD = clang

LD_LIBRARY_PATH := -L/usr/lib
INCL_PATH       := -Iinclude/

OUT_O_DIR := build
SRCS_DIRS := . methods
SRCS	    := $(wildcard $(addsuffix /*.c, $(SRCS_DIRS)))
OBJS	    := $(addprefix $(OUT_O_DIR)/,$(SRCS:.c=.o))
DEPS	    := $(OBJS:.o=.d)

ifdef DEBUG_INFO
	CFLAGS += -DDEBUG_INFO
endif

.PHONY: all clean
all: $(OUT_O_DIR)/main

$(OUT_O_DIR)/main: $(OBJS)
	$(LD) $^ -o $@ $(LDFLAGS) $(LD_LIBRARY_PATH)

$(OBJS): $(OUT_O_DIR)/%.o: %.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(INCL_PATH) -c $< -o $@

$(DEPS): $(OUT_O_DIR)/%.d: %.c
	@mkdir -p $(@D)
	$(CC) -E $(CFLAGS) $(INCL_PATH) $< -MM -MT $(@:.d=.o) > $@

clean:
	rm -rf $(OBJS) $(DEPS) $(OUT_O_DIR)

NODEPS = clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
include $(DEPS)
endif
