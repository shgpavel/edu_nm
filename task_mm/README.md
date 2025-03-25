## Overview
	Research task on vectorized matrix
	multiplication, comparing against MKL

## Status
	WIP

## Contains
	src/types/vector.c      -- updated vector lib
	src/types/matrix.c      -- updated matrix lib
	src/types/vector_avx.c  -- some extensions for vector
	src/methods/mult.c      -- actual matrix multiplication
	src/main.c              -- build matrix & test them

## Dependencies
	clang           >= 19.1.7
	jemalloc        >= 5.3
	OneAPI Basekit  >= 2025.0.1.46-1
	AVX CPU
	some libc libm

## Build
	source /opt/intel/oneapi/setvars.sh
	make

## Start
	./build/main
	python graph.py
