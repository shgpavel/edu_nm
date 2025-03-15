## Overview
	Research task on matrix multiplication
	comparing against MKL

## Status
	WIP

## Future plans
	1. FMA
	2. any matrix
	3. vtune
	4. better testing (graphics and random matrix generation)
	5. O(n^e) algo AlphaZero
	6. cblas interface

## Contains
	src/types/vector.c      -- updated vector lib
	src/types/matrix.c      -- updated matrix lib
	src/types/vector_avx.c  -- some extensions for vector
	src/methods/mm_ob.c     -- actual matrix multiplication
	src/main.c              -- build matrix & test them

## Dependencies
	clang           >= 19.1.7
	jemalloc        >= 5.3
	OneAPI Basekit  >= 2025.0.1.46-1
	AVX CPU
	some libc libm

## Build
	make
	
## Start
	source /opt/intel/oneapi/setvars.sh
	./build/main
	python graph.py
