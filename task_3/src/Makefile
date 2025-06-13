CFLAGS := -Wall -Wextra -O2 -std=c23 ${CFLAGS}

CC = clang
LD = clang

SRCS_DIRS = . lup matrix vector func
SRCS = $(wildcard $(addsuffix /*.c, $(SRCS_DIRS)))

LL = -lm -ljemalloc

main: $(SRCS)
	$(CC) $(SRCS) $(LL) $(CFLAGS)

clean:
	rm -rf a.out
