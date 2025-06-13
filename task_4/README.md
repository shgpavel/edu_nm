## Overview:
	If anyone wants to somehow use this code
	- Do not use this inv iter or use it only if you know what you doing
	- QR-method probably can be easly optimized by doing less allocations
		and my lecturer told me that I can't, but I still think tensor of
		Givens matrix in QR decomp can be optimized cause of it's
		eye-like matrix form
	- Iteration method is good

## Status:
    Released

## Contains:
	src/methods/gauss.c             -- just solving SLAE
	src/methods/inverse_iteration.c -- inv iter algo with two realizations

    reley shifts and regular from -norm_inf to +norm_inf 
	(the second one is "not that good")

    src/methods/power_method.c      -- simple power method
	src/methods/qr.c                -- QR method eigenvalues

    also contains hessenberg transform and givens rotations on hessenberg matrix

    src/methods/rng.h               -- simple time-based std "random"ng
    src/types/matrix.c              -- matrix methods
	src/types/vector.c              -- vector methods
	src/types/eigenpair.h           -- abstraction on vector and double
	src/common.h                    -- macros and defs
	src/main.c                      -- base logic

## Dependencies:
    jemalloc                        >= 5.3
    clang                           >= 15
    posix

## Build:
    make
    make clean
