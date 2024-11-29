## Overview:
	InterpolationSplines calculator, 100% educational
	Now fixed, all functions working (except ones mentioned in ISSUES.md)
	
## Status:
    Release

## Contains:
    src/methods/gauss.c             -- solving SLAE
    src/methods/polynoms.c          -- some polynoms operations like * +
	src/methods/lagrange.c          -- lagrange poly interpolation
    src/methods/newton.c            -- newton poly interpolation
	src/methods/linear_spline.c     -- linear spline interpolation (1 1)
    src/methods/quad_spline.c       -- spline 2 1
    src/methods/qube_spline.c       -- spline 3 2
    src/methods/penalty.c           -- spline 2 0, SLAE interp and newton2

    src/types/matrix.c              -- matrix methods
	src/types/vector.c              -- vector methods
    src/types/pair.h                -- simple pair with doubles

	(check drafter2 instead of this ->)
    src/draw/draw.c                 -- BAD cURL HTTP client
	src/draw/drawer.js              -- JS POST/GET server to draw with desmos

	src/funcs/funcs.h               -- functions from task, some aux
    src/common.h                    -- macros and defs
	src/main.c                      -- base logic

## Dependencies:
    jemalloc >= 5.3
    clang    >= 17
	libcurl  >= 8.7.1
	node     >= 21.7.3
	express  >= 4.19.2
	linux/gnu/xdg
	curl dependencies
    express dependencies
    libc math

## Build:
	cd src/draw
	npm install express
	cd ..
    make
    make clean

## Start:
	node src/draw/drawer.js
	./build/main
