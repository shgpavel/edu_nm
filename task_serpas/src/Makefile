CFLAGS := -Wall -Wextra -std=c23 -fopenmp -flto $(CFLAGS)
LDFLAGS := -lm -lomp $(LDFLAGS)

-fsanitize=address -fsanitize=undefined

CFLAGS += -O4 -ffast -ffast-mem -ffast-math -ffp-contract=fast -fexcess-precision=fast

CFLAGS += -fswp-maxopers=2000 -fcache-opt -funroll-loops -flist-prefetch

#CFLAGS += -O3 -mavx -mfma -funroll-loops

CC = lcc
#CC = clang

LD_LIBRARY_PATH = -L/usr/lib
INCL_PATH = -Iinclude/

all: sca25lib.so

sca25lib.so: sca25lib.c
	$(CC) -fPIC -shared $(CFLAGS) $(LDFLAGS) $(INCL_PATH) $< -o $@

test: test.c
	$(CC) $(CFLAGS) test.c

clean:
	rm -f sca25lib.o sca25lib.so a.out
