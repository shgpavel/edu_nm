CFLAGS := -Wall -Wextra -Wconversion -Wformat -Wformat-security \
-O2 -std=c23 ${CFLAGS}
# -fsanitize=address -fsanitize=undefined 

LDFLAGS := -lm -ljemalloc -lcurl ${LDFLAGS}
# -fsanitize=address -fsanitize=undefined

CC = clang
LD = clang

OUT_O_DIR = build
SRCS_DIRS = . methods types draw
SRCS = $(wildcard $(addsuffix /*.c, $(SRCS_DIRS)))
OBJS := $(addprefix $(OUT_O_DIR)/,$(SRCS:.c=.o))
DEPS = $(OBJS:.o=.d)

.PHONY: all
all: $(OUT_O_DIR)/main

$(OUT_O_DIR)/main: $(OBJS)
	$(LD) $^ -o $@ $(LDFLAGS)

$(OBJS): $(OUT_O_DIR)/%.o: %.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(DEPS): $(OUT_O_DIR)/%.d: %.c
	@mkdir -p $(@D)
	$(CC) -E $(CFLAGS) $< -MM -MT $(@:.d=.o) > $@

.PHONY: clean
clean:
	rm -rf $(OBJS) $(DEPS) $(OUT_O_DIR)

NODEPS = clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
include $(DEPS)
endif
