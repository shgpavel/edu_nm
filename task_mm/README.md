## Overview
	Research task on matrix multiplication
	Comparing with OpenBLAS

## Status
	Started

## Contains
	src/types/vector.c      -- updated vector lib
	src/types/matrix.c      -- updated matrix lib
	src/types/vector_avx.c  -- some extensions for vector
	src/main.c              -- build matrix & test

## Dependencies
	clang     >= 19.1.7
	jemalloc  >= 5.3
	OpenBLAS  >= 0.3.28
	Intel CPU >= Haswell
	libc libm

## Build
	make
	
## Start
	./build/main
