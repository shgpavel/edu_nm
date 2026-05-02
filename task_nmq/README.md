## Overview

This directory includes my submissions for the **NM Quad NPU** programming
tournament. Targets the K1879VM8Ya SoC (NMC4 VLIW DSP cores) on the NM Quad
PCIe accelerator.

Preferred language: NMC assembly (`.S`). C is used only inside harnesses.

## Status
	Release

## Tasks
	src/a+b/      -- 32-bit signed addition
	src/demomm/   -- int32 matrix multiplication
	src/fastpow/  -- modular exponentiation
	src/fft128/   -- 128-point FFT
	src/fft512/   -- 512-point FFT
	src/lastfct/  -- last non-zero digit of n!
	src/rc4/      -- RC4 stream cipher
	src/slae/     -- system of linear equations (Gaussian elimination)

	desc/<task>/main.txt -- problem statements (in Russian)

## Dependencies
	NMC-SDK 1.4.249 (provides nmc-gcc, nmc-qemu, nmc-objdump, headers, libs)
	GNU Make
	libc, shell

## Configuration

`NMC_PREFIX` points at the SDK install root (the directory that contains
`bin/nmc-gcc`, `nmc/include/`, `nmc/lib/nmc4/`). Each task's Makefile
references it; nothing else carries an absolute path.

Set it in your shell before building:

```sh
export NMC_PREFIX=/path/to/NMC-SDK-1.4.249.d4227c9.x86_64/usr/local
```

Or pass it on the command line:

```sh
make NMC_PREFIX=/path/to/NMC-SDK-.../usr/local
```

The Makefile default is a placeholder and will not build until you
override it.

## Build
	cd src/<task>
	make

## Test
	make test

Runs the open test vectors through `nmc-qemu`.

## Bench
	make bench

Wall-clock host timing across the QEMU run. Relative proxy only,
since QEMU does not model NMC4 cycle timing.
