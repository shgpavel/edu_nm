## Overview
	Research task on matrix multiplication
	Comparing with OpenBLAS

## Status
	In progress

## Contains
	src/types/vector.c      -- updated vector lib
	src/types/matrix.c      -- updated matrix lib
	src/types/vector_avx.c  -- some extensions for vector
	src/methods/mm_ob.c     -- actual matrix multiplication
	src/main.c              -- build matrix & test them

## Dependencies
	clang     >= 19.1.7
	jemalloc  >= 5.3
	OpenBLAS  >= 0.3.28
	AVX CPU
	libc libm

## Build
	make
	
## Start
	taskset -c 0 ./build/main
	python graph.py
