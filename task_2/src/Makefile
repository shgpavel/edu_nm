CFLAGS := -Wall -Wextra -O2 -std=c++23 ${CFLAGS}

CC = clang++
LD = clang++

SRCS_DIRS = . MPI MS LUP QR EigenWrap Funcs
SRCS = $(wildcard $(addsuffix /*.cpp, $(SRCS_DIRS)))

main: main.cpp
	$(CC) $(SRCS) $(CFLAGS)

clean:
	rm -rf a.out 
