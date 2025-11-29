## overview

This directory includes my submissions for a September 2025 programming competition
centered on the Elbrus-2000 (e2k) CPU architecture.
The contest brought together 31 participants and focused on practical
linear-algebra workloads: matrix multiplication, determinant calculation, and solving systems
of equations. Participants used C/C++ (asm, intrinsics as well) and worked through an SSH-based setup that provided
access to a real Elbrus machine and a cross-compiler targeting e2k.

https://habr.com/ru/articles/959742/

https://www.serpas.ru/contest/sca25/results/qualification.html

https://www.serpas.ru/contest/sca25/results/task1.html

https://www.serpas.ru/contest/sca25/results/task2.html

https://www.serpas.ru/contest/sca25/results/task3.html

https://www.youtube.com/live/xqBdo9USor8

https://youtu.be/i5uQbg1xJO4

## status
        release

## contains
        src/sca25lib.c -- all methods
        src/test.c     -- my tests
        src/stest_x86  -- championship's scoring system (x86 variant)

## common required
        libc
        libm
        GNU Make
        OpenMP

### e2k
        Elbrus Linux
        lcc      >= 1.29.12

### x86_64
        clang    >= 21.1.6

## build
        make
builds shared object for stest or stest_x86

## start
        ./stest_x86 -t[1,2,3] -r[int]
