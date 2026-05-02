	.file	"solution.c"
	.text
	.align 1;
	.macro FFT128_BUTTERFLY
	ar3 = ar0;
	ar4 = ar1;
	ar5 = ar2;
	fpu 0 rep vlen vreg0 = [ar4++];
	fpu 0 rep vlen vreg2 = [ar5++];
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar4++];
	fpu 0 rep vlen vreg3 = [ar5++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 .double vreg4 = vreg0 * vreg2;
	fpu 0 .double vreg5 = vreg1 * vreg2;
	fpu 0 .double vreg4 = vreg1 * vreg3 + vreg4;
	fpu 0 .double vreg5 = -vreg0 * vreg3 + vreg5;
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	ar3 = ar0;
	ar4 = ar1;
	fpu 0 rep vlen [ar3++] = vreg0;
	fpu 0 rep vlen [ar3++] = vreg1;
	fpu 0 rep vlen [ar4++] = vreg2;
	fpu 0 rep vlen [ar4++] = vreg3;
	.endm

	.macro FFT128_BUTTERFLY_PRELOAD
	ar3 = ar0;
	ar4 = ar1;
	fpu 0 rep vlen vreg0 = [ar4++];
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar4++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 .double vreg2 = vreg0 * vreg4;
	fpu 0 .double vreg3 = vreg1 * vreg4;
	fpu 0 .double vreg2 = vreg1 * vreg5 + vreg2;
	fpu 0 .double vreg3 = -vreg0 * vreg5 + vreg3;
	fpu 0 .double vreg0 = vreg6 + vreg2;
	fpu 0 .double vreg1 = vreg7 + vreg3;
	fpu 0 .double vreg2 = vreg6 - vreg2;
	fpu 0 .double vreg3 = vreg7 - vreg3;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	.endm

	.macro FFT128X16_BUTTERFLY
	ar3 = ar0;
	ar4 = ar1;
	ar5 = ar2;
	fpu 0 rep vlen vreg0 = [ar4++];
	fpu 0 rep vlen vreg1 = [ar4++];
	fpu 0 rep vlen vreg2 = [ar5];
	ar5 = ar5 + 2;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg3 = [ar5];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 .double vreg4 = vreg0 * vreg2;
	fpu 0 .double vreg5 = vreg1 * vreg2;
	fpu 0 .double vreg4 = vreg1 * vreg3 + vreg4;
	fpu 0 .double vreg5 = -vreg0 * vreg3 + vreg5;
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	.endm

	.macro FFT128X16_BUTTERFLY_PRELOAD
	ar3 = ar0;
	ar4 = ar1;
	fpu 0 rep vlen vreg0 = [ar4++];
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar4++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 .double vreg2 = vreg0 * vreg4;
	fpu 0 .double vreg3 = vreg1 * vreg4;
	fpu 0 .double vreg2 = vreg1 * vreg5 + vreg2;
	fpu 0 .double vreg3 = -vreg0 * vreg5 + vreg3;
	fpu 0 .double vreg0 = vreg6 + vreg2;
	fpu 0 .double vreg1 = vreg7 + vreg3;
	fpu 0 .double vreg2 = vreg6 - vreg2;
	fpu 0 .double vreg3 = vreg7 - vreg3;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	.endm

	.macro FFT128_ROWX16_GATHER_LANE
	ar4,gr4 = [ar5++];
	[ar1++] = ar4;
	[ar2++] = gr4;
	ar5 = ar5 + -258;
	.endm

	.macro FFT128_ROWX16_GATHER_LAST
	ar4,gr4 = [ar5++];
	[ar1++] = ar4;
	[ar2++] = gr4;
	.endm

	.macro FFT128_ROWX16_SCATTER_LANE
	ar4,gr4 = [ar0++];
	[ar5++] = ar4,gr4;
	ar4,gr4 = [ar2++];
	[ar5++] = ar4,gr4;
	ar5 = ar5 + 508;
	.endm

	.macro FFT128_ROWX16_SCATTER_LAST
	ar4,gr4 = [ar0++];
	[ar5++] = ar4,gr4;
	ar4,gr4 = [ar2++];
	[ar5++] = ar4,gr4;
	.endm

	.if 0
	.globl _fft1d_128
	.type	_fft1d_128, @function
_fft1d_128:
	gr7 = [ar7 + -3];
	push ar6,gr6;
	push ar5,gr5;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	ar7 = ar7 + 2;
	[ar7 + -1] = gr7;
	ar0 = gr7;
	gr5 = ar0;
	ar2 = _g_bitrev;
	ar1 = _g_fft_tmp;
	gr6 = 0;
	vlen = gr6;
	gr1 = 32;
	.align;
LFFT128_BITREV_STAGE12:
	gr0 = [ar2++];
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg0 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar3++];
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg2 = [ar3++];
	fpu 0 rep vlen vreg3 = [ar3++];
	fpu 0 .double vreg4 = vreg0 + vreg2;
	fpu 0 .double vreg5 = vreg1 + vreg3;
	fpu 0 .double vreg6 = vreg0 - vreg2;
	fpu 0 .double vreg7 = vreg1 - vreg3;
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg0 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar3++];
	gr0 = gr5 + gr0;
	ar4 = gr0;
	ar3 = ar1;
	fpu 0 rep vlen vreg2 = [ar4++];
	fpu 0 rep vlen vreg3 = [ar4++];
	fpu 0 .double vreg0 = vreg0 + vreg2;
	fpu 0 .double vreg1 = vreg1 + vreg3;
	fpu 0 .double vreg2 = vreg2 + vreg2;
	fpu 0 .double vreg3 = vreg3 + vreg3;
	fpu 0 .double vreg2 = vreg0 - vreg2;
	fpu 0 .double vreg3 = vreg1 - vreg3;
	fpu 0 .double vreg0 = vreg4 + vreg0;
	fpu 0 .double vreg1 = vreg5 + vreg1;
	fpu 0 .double vreg4 = vreg4 + vreg4;
	fpu 0 .double vreg5 = vreg5 + vreg5;
	fpu 0 .double vreg4 = vreg4 - vreg0;
	fpu 0 .double vreg5 = vreg5 - vreg1;
	ar4 = ar1 + 8;
	fpu 0 rep vlen [ar3++] = vreg0;
	fpu 0 rep vlen [ar3++] = vreg1;
	fpu 0 rep vlen [ar4++] = vreg4;
	fpu 0 rep vlen [ar4++] = vreg5;
	fpu 0 .double vreg3 = vreg6 + vreg3;
	fpu 0 .double vreg2 = vreg7 - vreg2;
	fpu 0 .double vreg6 = vreg6 + vreg6;
	fpu 0 .double vreg7 = vreg7 + vreg7;
	fpu 0 .double vreg6 = vreg6 - vreg3;
	fpu 0 .double vreg7 = vreg7 - vreg2;
	ar3 = ar1 + 4;
	ar4 = ar1 + 12;
	fpu 0 rep vlen [ar3++] = vreg3;
	fpu 0 rep vlen [ar3++] = vreg2;
	fpu 0 rep vlen [ar4++] = vreg6;
	fpu 0 rep vlen [ar4++] = vreg7;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_BITREV_STAGE12;
	ar1 = ar1 + 16;
	nul;
	set fp_wait;
	set fp_branch;
	ar0 = _g_fft_tmp;
	gr1 = 16;
	.align;
LFFT128_STAGE3_R0:
	ar3 = ar0;
	ar4 = ar0 + 16;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE3_R0;
	ar0 = ar0 + 32;
	nul;
	ar0 = _g_fft_tmp + 4;
	ar1 = _g_fft_tmp + 20;
	ar2 = _g_roots + 64;
	gr1 = 16;
	fpu 0 rep vlen vreg4 = [ar2++];
	fpu 0 rep vlen vreg5 = [ar2++];
	.align;
LFFT128_STAGE3_R16:
	FFT128_BUTTERFLY_PRELOAD
	ar1 = ar1 + 32;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE3_R16;
	ar0 = ar0 + 32;
	nul;
	ar0 = _g_fft_tmp + 8;
	gr1 = 16;
	.align;
LFFT128_STAGE3_R32:
	ar3 = ar0;
	ar4 = ar0 + 16;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg5;
	fpu 0 .double vreg1 = vreg7 - vreg4;
	fpu 0 .double vreg2 = vreg6 - vreg5;
	fpu 0 .double vreg3 = vreg7 + vreg4;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE3_R32;
	ar0 = ar0 + 32;
	nul;
	ar0 = _g_fft_tmp + 12;
	ar1 = _g_fft_tmp + 28;
	ar2 = _g_roots + 192;
	gr1 = 16;
	fpu 0 rep vlen vreg4 = [ar2++];
	fpu 0 rep vlen vreg5 = [ar2++];
	.align;
LFFT128_STAGE3_R48:
	FFT128_BUTTERFLY_PRELOAD
	ar1 = ar1 + 32;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE3_R48;
	ar0 = ar0 + 32;
	nul;
	set fp_wait;
	set fp_branch;
	ar0 = _g_fft_tmp;
	gr1 = 8;
	.align;
LFFT128_STAGE4_R0:
	ar3 = ar0;
	ar4 = ar0 + 32;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE4_R0;
	ar0 = ar0 + 64;
	nul;
	ar0 = _g_fft_tmp + 4;
	ar1 = _g_fft_tmp + 36;
	ar2 = _g_roots + 32;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2++];
	fpu 0 rep vlen vreg5 = [ar2++];
	.align;
LFFT128_STAGE4_R8:
	FFT128_BUTTERFLY_PRELOAD
	ar1 = ar1 + 64;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE4_R8;
	ar0 = ar0 + 64;
	nul;
	ar0 = _g_fft_tmp + 8;
	ar1 = _g_fft_tmp + 40;
	ar2 = _g_roots + 64;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2++];
	fpu 0 rep vlen vreg5 = [ar2++];
	.align;
LFFT128_STAGE4_R16:
	FFT128_BUTTERFLY_PRELOAD
	ar1 = ar1 + 64;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE4_R16;
	ar0 = ar0 + 64;
	nul;
	ar0 = _g_fft_tmp + 12;
	ar1 = _g_fft_tmp + 44;
	ar2 = _g_roots + 96;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2++];
	fpu 0 rep vlen vreg5 = [ar2++];
	.align;
LFFT128_STAGE4_R24:
	FFT128_BUTTERFLY_PRELOAD
	ar1 = ar1 + 64;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE4_R24;
	ar0 = ar0 + 64;
	nul;
	ar0 = _g_fft_tmp + 16;
	gr1 = 8;
	.align;
LFFT128_STAGE4_R32:
	ar3 = ar0;
	ar4 = ar0 + 32;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg5;
	fpu 0 .double vreg1 = vreg7 - vreg4;
	fpu 0 .double vreg2 = vreg6 - vreg5;
	fpu 0 .double vreg3 = vreg7 + vreg4;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE4_R32;
	ar0 = ar0 + 64;
	nul;
	ar0 = _g_fft_tmp + 20;
	ar1 = _g_fft_tmp + 52;
	ar2 = _g_roots + 160;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2++];
	fpu 0 rep vlen vreg5 = [ar2++];
	.align;
LFFT128_STAGE4_R40:
	FFT128_BUTTERFLY_PRELOAD
	ar1 = ar1 + 64;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE4_R40;
	ar0 = ar0 + 64;
	nul;
	ar0 = _g_fft_tmp + 24;
	ar1 = _g_fft_tmp + 56;
	ar2 = _g_roots + 192;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2++];
	fpu 0 rep vlen vreg5 = [ar2++];
	.align;
LFFT128_STAGE4_R48:
	FFT128_BUTTERFLY_PRELOAD
	ar1 = ar1 + 64;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE4_R48;
	ar0 = ar0 + 64;
	nul;
	ar0 = _g_fft_tmp + 28;
	ar1 = _g_fft_tmp + 60;
	ar2 = _g_roots + 224;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2++];
	fpu 0 rep vlen vreg5 = [ar2++];
	.align;
LFFT128_STAGE4_R56:
	FFT128_BUTTERFLY_PRELOAD
	ar1 = ar1 + 64;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_STAGE4_R56;
	ar0 = ar0 + 64;
	nul;
	set fp_wait;
	gr0 = 64;
	gr2 = 16;
	gr4 = 16;
	gr7 = 4;
	set fp_branch;
	.align;
LFFT128_STAGE:
	ar0 = _g_fft_tmp;
	gr3 = gr7;
	.align;
LFFT128_BLOCK_SETUP:
	ar2 = _g_roots;
	ar1 = ar0 + gr0;
	gr1 = gr4;
LFFT128_INNER:
	FFT128_BUTTERFLY
	ar2 = ar2 + gr2 with gr1 = gr1 - 1;
	ar1 = ar1 + 4;
	if <>0 delayed goto LFFT128_INNER;
	ar0 = ar0 + 4;
	nul;
	gr3 = gr3 - 1;
	if =0 delayed goto LFFT128_STAGE_DONE;
	nul;
	nul;
	ar0 = ar1;
	delayed goto LFFT128_BLOCK_SETUP;
	nul;
	nul;
	.align;
LFFT128_STAGE_DONE:
	set fp_wait;
	gr0 = gr0 + gr0;
	gr2 = gr2 >> 1;
	gr4 = gr4 + gr4;
	gr7 = gr7 >> 1;
	gr7;
	if =0 delayed goto LFFT128_COPY_BACK_SETUP;
	nul;
	nul;
	set fp_branch;
	delayed goto LFFT128_STAGE;
	nul;
	nul;
LFFT128_COPY_BACK_SETUP:
	set fp_wait;
	gr6 = [ar7 + -1];
	gr0 = _g_column;
	gr6 - gr0;
	if =0 delayed goto LFFT128_RET;
	nul;
	nul;
	ar0 = gr6;
	ar1 = _g_fft_tmp;
	gr1 = 128;
LFFT128_COPY_BACK:
	ar5,gr5 = [ar1++];
	[ar0++] = ar5,gr5;
	ar5,gr5 = [ar1++];
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128_COPY_BACK;
	[ar0++] = ar5,gr5;
	nul;
LFFT128_RET:
	ar7 = ar7 + -2;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size	_fft1d_128, .-_fft1d_128
	.endif

	.align 1;
	.globl _fft1d_128x16
	.type	_fft1d_128x16, @function
_fft1d_128x16:
	push ar6,gr6;
	push ar5,gr5;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	gr6 = 15;
	vlen = gr6;
	set fp_branch;
	ar0 = _g_col_batch;
	ar2 = _g_bitrev;
	ar1 = _g_col_tmp;
	gr5 = ar0;
	gr1 = 32;
	.align;
LFFT128X16_BITREV_STAGE12:
	gr0 = [ar2++];
	gr0 = gr0 << 4;
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg0 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar3++];
	gr0 = gr0 << 4;
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg2 = [ar3++];
	fpu 0 rep vlen vreg3 = [ar3++];
	fpu 0 .double vreg4 = vreg0 + vreg2;
	fpu 0 .double vreg5 = vreg1 + vreg3;
	fpu 0 .double vreg6 = vreg0 - vreg2;
	fpu 0 .double vreg7 = vreg1 - vreg3;
	gr0 = gr0 << 4;
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg0 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar3++];
	gr0 = gr0 << 4;
	gr0 = gr5 + gr0;
	ar4 = gr0;
	fpu 0 rep vlen vreg2 = [ar4++];
	fpu 0 rep vlen vreg3 = [ar4++];
	fpu 0 .double vreg0 = vreg0 + vreg2;
	fpu 0 .double vreg1 = vreg1 + vreg3;
	fpu 0 .double vreg2 = vreg2 + vreg2;
	fpu 0 .double vreg3 = vreg3 + vreg3;
	fpu 0 .double vreg2 = vreg0 - vreg2;
	fpu 0 .double vreg3 = vreg1 - vreg3;
	fpu 0 .double vreg0 = vreg4 + vreg0;
	fpu 0 .double vreg1 = vreg5 + vreg1;
	fpu 0 .double vreg4 = vreg4 + vreg4;
	fpu 0 .double vreg5 = vreg5 + vreg5;
	fpu 0 .double vreg4 = vreg4 - vreg0;
	fpu 0 .double vreg5 = vreg5 - vreg1;
	ar3 = ar1;
	ar4 = ar1 + 128;
	fpu 0 rep vlen [ar1++] = vreg0;
	fpu 0 rep vlen [ar1++] = vreg1;
	fpu 0 rep vlen [ar4++] = vreg4;
	fpu 0 rep vlen [ar4++] = vreg5;
	fpu 0 .double vreg3 = vreg6 + vreg3;
	fpu 0 .double vreg2 = vreg7 - vreg2;
	fpu 0 .double vreg6 = vreg6 + vreg6;
	fpu 0 .double vreg7 = vreg7 + vreg7;
	fpu 0 .double vreg6 = vreg6 - vreg3;
	fpu 0 .double vreg7 = vreg7 - vreg2;
	ar3 = ar3 + 64;
	ar4 = ar3 + 128;
	fpu 0 rep vlen [ar3++] = vreg3;
	fpu 0 rep vlen [ar3++] = vreg2;
	fpu 0 rep vlen [ar4++] = vreg6;
	fpu 0 rep vlen [ar4++] = vreg7;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_BITREV_STAGE12;
	ar1 = ar1 + 192;
	nul;
	set fp_wait;
	set fp_branch;
	ar0 = _g_col_tmp;
	gr1 = 16;
	.align;
LFFT128X16_STAGE3_R0:
	ar3 = ar0;
	ar4 = ar0 + 256;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE3_R0;
	ar0 = ar0 + 512;
	nul;
	set fp_wait;
	ar2 = _g_roots + 64;
	gr1 = 16;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 64;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 320;
	set fp_branch;
	.align;
LFFT128X16_STAGE3_R16:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 512;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE3_R16;
	ar0 = ar0 + 512;
	nul;
	set fp_wait;
	ar0 = _g_col_tmp + 128;
	gr1 = 16;
	set fp_branch;
	.align;
LFFT128X16_STAGE3_R32:
	ar3 = ar0;
	ar4 = ar0 + 256;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg5;
	fpu 0 .double vreg1 = vreg7 - vreg4;
	fpu 0 .double vreg2 = vreg6 - vreg5;
	fpu 0 .double vreg3 = vreg7 + vreg4;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE3_R32;
	ar0 = ar0 + 512;
	nul;
	set fp_wait;
	ar2 = _g_roots + 192;
	gr1 = 16;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 192;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 448;
	set fp_branch;
	.align;
LFFT128X16_STAGE3_R48:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 512;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE3_R48;
	ar0 = ar0 + 512;
	nul;
	set fp_wait;
	set fp_branch;
	ar0 = _g_col_tmp;
	gr1 = 8;
	.align;
LFFT128X16_STAGE4_R0:
	ar3 = ar0;
	ar4 = ar0 + 512;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE4_R0;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 32;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 64;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 576;
	set fp_branch;
	.align;
LFFT128X16_STAGE4_R8:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE4_R8;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 64;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 128;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 640;
	set fp_branch;
	.align;
LFFT128X16_STAGE4_R16:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE4_R16;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 96;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 192;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 704;
	set fp_branch;
	.align;
LFFT128X16_STAGE4_R24:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE4_R24;
	ar0 = ar0 + 1024;
	nul;
	ar0 = _g_col_tmp + 256;
	gr1 = 8;
	set fp_branch;
	.align;
LFFT128X16_STAGE4_R32:
	ar3 = ar0;
	ar4 = ar0 + 512;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg5;
	fpu 0 .double vreg1 = vreg7 - vreg4;
	fpu 0 .double vreg2 = vreg6 - vreg5;
	fpu 0 .double vreg3 = vreg7 + vreg4;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE4_R32;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 160;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 320;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 832;
	set fp_branch;
	.align;
LFFT128X16_STAGE4_R40:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE4_R40;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 192;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 384;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 896;
	set fp_branch;
	.align;
LFFT128X16_STAGE4_R48:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE4_R48;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 224;
	gr1 = 8;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 448;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 960;
	set fp_branch;
	.align;
LFFT128X16_STAGE4_R56:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE4_R56;
	ar0 = ar0 + 1024;
	nul;
	set fp_wait;
	set fp_branch;
	ar0 = _g_col_tmp;
	gr1 = 4;
	.align;
LFFT128X16_STAGE5_R0:
	ar3 = ar0;
	ar4 = ar0 + 1024;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R0;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 16;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 64;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1088;
	.align;
LFFT128X16_STAGE5_R4:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R4;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 32;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 128;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1152;
	.align;
LFFT128X16_STAGE5_R8:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R8;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 48;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 192;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1216;
	.align;
LFFT128X16_STAGE5_R12:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R12;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 64;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 256;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1280;
	.align;
LFFT128X16_STAGE5_R16:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R16;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 80;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 320;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1344;
	.align;
LFFT128X16_STAGE5_R20:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R20;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 96;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 384;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1408;
	.align;
LFFT128X16_STAGE5_R24:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R24;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 112;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 448;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1472;
	.align;
LFFT128X16_STAGE5_R28:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R28;
	ar0 = ar0 + 2048;
	nul;
	ar0 = _g_col_tmp + 512;
	gr1 = 4;
	.align;
LFFT128X16_STAGE5_R32:
	ar3 = ar0;
	ar4 = ar0 + 1024;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg5;
	fpu 0 .double vreg1 = vreg7 - vreg4;
	fpu 0 .double vreg2 = vreg6 - vreg5;
	fpu 0 .double vreg3 = vreg7 + vreg4;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R32;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 144;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 576;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1600;
	.align;
LFFT128X16_STAGE5_R36:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R36;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 160;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 640;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1664;
	.align;
LFFT128X16_STAGE5_R40:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R40;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 176;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 704;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1728;
	.align;
LFFT128X16_STAGE5_R44:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R44;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 192;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 768;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1792;
	.align;
LFFT128X16_STAGE5_R48:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R48;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 208;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 832;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1856;
	.align;
LFFT128X16_STAGE5_R52:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R52;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 224;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 896;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1920;
	.align;
LFFT128X16_STAGE5_R56:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R56;
	ar0 = ar0 + 2048;
	nul;
	ar2 = _g_roots + 240;
	gr1 = 4;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 960;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 1984;
	.align;
LFFT128X16_STAGE5_R60:
	FFT128X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 2048;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_STAGE5_R60;
	ar0 = ar0 + 2048;
	nul;
	set fp_wait;
	gr0 = 2048;
	gr2 = 8;
	gr4 = 32;
	gr7 = 2;
	set fp_branch;
	.align;
LFFT128X16_STAGE:
	ar0 = _g_col_tmp;
	gr3 = gr7;
	.align;
LFFT128X16_BLOCK_SETUP:
	ar2 = _g_roots;
	ar1 = ar0 + gr0;
	gr1 = gr4;
LFFT128X16_INNER:
	FFT128X16_BUTTERFLY
	ar2 = ar2 + gr2 with gr1 = gr1 - 1;
	ar1 = ar1 + 64;
	if <>0 delayed goto LFFT128X16_INNER;
	ar0 = ar0 + 64;
	nul;
	gr3 = gr3 - 1;
	if =0 delayed goto LFFT128X16_STAGE_DONE;
	nul;
	nul;
	ar0 = ar1;
	delayed goto LFFT128X16_BLOCK_SETUP;
	nul;
	nul;
	.align;
LFFT128X16_STAGE_DONE:
	set fp_wait;
	gr0 = gr0 + gr0;
	gr2 = gr2 >> 1;
	gr4 = gr4 + gr4;
	gr7 = gr7 >> 1;
	gr7;
	if =0 delayed goto LFFT128X16_RET;
	gr6 = gr7 + -1;
	set fp_branch;
	if <>0 delayed goto LFFT128X16_STAGE;
	ar5 = [_g_x16_col_out];
	ar0 = _g_col_tmp;
	gr6 = ar5;
	gr6;
	if =0 delayed goto LFFT128X16_STAGE;
	ar1 = ar0 + gr0;
	ar2 = _g_roots;
LFFT128X16_COL_FINAL:
	ar6 = ar5 + 16384;
	gr1 = 64;
LFFT128X16_COL_FINAL_INNER:
	fpu 0 rep vlen vreg0 = [ar1++];
	fpu 0 rep vlen vreg1 = [ar1++];
	fpu 0 rep vlen vreg2 = [ar2];
	ar2 = ar2 + 2;
	fpu 0 rep vlen vreg6 = [ar0++];
	fpu 0 rep vlen vreg3 = [ar2];
	ar2 = ar2 + 2;
	fpu 0 rep vlen vreg7 = [ar0++];
	fpu 0 .double vreg4 = vreg0 * vreg2;
	fpu 0 .double vreg5 = vreg1 * vreg2;
	fpu 0 .double vreg4 = vreg1 * vreg3 + vreg4;
	fpu 0 .double vreg5 = -vreg0 * vreg3 + vreg5;
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	fpu 0 .packer = (vreg0,vreg1) with .float <= .double;
	fpu rep 16 [ar5++] = .packer;
	ar5 = ar5 + 224;
	fpu 0 .packer = (vreg2,vreg3) with .float <= .double;
	fpu rep 16 [ar6++] = .packer;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT128X16_COL_FINAL_INNER;
	ar6 = ar6 + 224;
	nul;
	[_g_x16_col_out] = gr1;
	set fp_wait;
	nul;
LFFT128X16_RET:
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size	_fft1d_128x16, .-_fft1d_128x16
	.align 1;
.globl _contest_fft_2d
	.type	_contest_fft_2d, @function
_contest_fft_2d:
	push ar6,gr6;
	push ar5,gr5;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	ar7 = ar7 + 64;
	ar2 = ar7 + -80;
	gr7 = [ar2 + -1];
	[ar7 + -40] = gr7;
	ar0 = gr7;
	gr1 = 1023;
	gr0 = 0;
	gr2 = 0;
LZERO_SCAN128:
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	gr1 = gr1 + -1;
	if <>0 delayed goto LZERO_SCAN128_PREF;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
LZERO_SCAN128_PREF:
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	gr1 = gr1 + -1;
	if <>0 delayed goto LZERO_SCAN128_PREF;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	ar5,gr5 = [ar0++] with gr2 = gr2 or gr0;
	gr0 = ar5;
	gr0 = gr0 or gr5;
	gr2 = gr2 or gr0;
	gr2;
	if <>0 goto LZERO_NONZERO128;
	goto LRET_CONTEST;
LZERO_NONZERO128:
	gr6 = [_g_tables_ready];
	gr6;
	if <>0 delayed goto L19;
	nul;
	nul;
	ar3 = _g_roots;
	ar2 = gr6;
	ar0 = gr6;
	ar1 = gr6;
	gr1 = 1072693248;
	[ar7 + -20] = gr6;
	[ar7 + -19] = gr1;
	gr6 = ar3;
	[ar7 + -28] = gr2;
	[ar7 + -27] = gr3;
	[ar7 + -36] = gr4;
	[ar7 + -35] = gr5;
	gr5 = ar2;
	ar6 = ar7 + -12;
	gr1 = [ar6];
	[ar7 + -48] = gr1;
	ar4 = ar7 + -11;
	gr1 = [ar4];
	[ar7 + -47] = gr1;
	ar5 = [ar7 + -20];
	ar4 = gr6;
	[ar4] = ar5;
	ar5 = [ar7 + -19];
	[ar4 + 1] = ar5;
	ar6 = ar4 + 2;
	[ar6] = ar0;
	[ar4 + 3] = ar1;
	gr3 = gr5;
	ar6 = 7;
	gr2 = 0;
	gr4 = 1;
L23:
	gr2 = gr2 + gr2;
	gr7 = gr3 and gr4;
	gr2 = gr7 or gr2;
	gr3 = gr3 >> 1;
	ar6 = ar6 + -1;
	gr0 = ar6;
	gr0;
	if <>0 goto L23;
	gr1 = _g_bitrev;
	gr3 = gr1 + gr5;
	ar4 = gr3;
	gr2 = gr2 << 2;
	[ar4] = gr2;
	ar5 = gr5;
	ar5 = ar5 + 1;
	gr5 = ar5;
	gr2 = 128;
	gr5 - gr2;
	if =0 goto L36;
	gr2 = 63;
	gr3 = gr5 and gr2;
	gr3 = gr3 + gr3;
	gr7 = _k_base_cos;
	gr7 = gr7 + gr3;
	ar6 = gr7;
	gr1 = [ar6];
	[ar7 + -20] = gr1;
	gr1 = [ar6 + 1];
	[ar7 + -19] = gr1;
	gr7 = _k_base_sin;
	gr3 = gr7 + gr3;
	ar6 = gr3;
	ar0 = [ar6];
	ar1 = [ar6 + 1];
	gr2 - gr5;
	if u>= goto L21;
	ar6 = gr1;
	ar3 = ar6 + -2147483648;
	[ar7 + -19] = ar3;
	[ar7 + -56] = ar0;
	ar6 = ar1 + -2147483648;
	[ar7 + -55] = ar6;
	ar4 = ar7 + -56;
	ar0 = [ar4];
	ar1 = ar6;
L21:
	ar6 = gr6;
	ar6 = ar6 + 4;
	gr6 = ar6;
	ar5 = [ar7 + -20];
	ar4 = gr6;
	[ar4] = ar5;
	ar5 = [ar7 + -19];
	[ar4 + 1] = ar5;
	ar6 = ar4 + 2;
	[ar6] = ar0;
	[ar4 + 3] = ar1;
	gr3 = gr5;
	ar6 = 7;
	gr2 = 0;
	gr4 = 1;
	goto L23;
L36:
	ar6 = ar7 + -48;
	gr1 = [ar6];
	[ar7 + -12] = gr1;
	ar3 = ar7 + -47;
	gr1 = [ar3];
	[ar7 + -11] = gr1;
	gr6 = 1;
	[_g_tables_ready] = gr6;
L19:
	set fp_branch;
	ar2 = ar7 + -40;
	ar3 = [ar2];
	ar6 = _g_row_scratch;
	ar2 = ar7 + -1;
	[ar7 + -28] = ar2;
	ar4 = ar7 + -12;
	gr1 = [ar4];
	[ar7 + -36] = gr1;
	ar4 = ar7 + -11;
	gr1 = [ar4];
	[ar7 + -35] = gr1;
	gr5 = 8;
L28:
	ar3 = ar3 + 3840;
	ar0 = _g_col_batch;
	gr6 = 128;
LROWX16_GATHER_POINT:
	ar5 = ar3;
	ar1 = _g_column_f;
	ar2 = _g_column_f + 16;
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LANE
	FFT128_ROWX16_GATHER_LAST
	ar1 = _g_column_f;
	ar2 = ar0;
	fpu rep 8 .packer = [ar1++] with .double <= .float;
	fpu rep 16 [ar2++] = .packer;
	ar1 = _g_column_f + 16;
	fpu rep 8 .packer = [ar1++] with .double <= .float;
	fpu rep 16 [ar2++] = .packer;
	ar0 = ar2;
	gr6 = gr6 + -1;
	gr6;
	if <>0 delayed goto LROWX16_GATHER_POINT;
	ar3 = ar3 + 2;
	nul;
	set fp_wait;
	call _fft1d_128x16;
	ar0 = _g_col_tmp;
	ar1 = ar6;
	gr6 = 128;
LROWX16_SCATTER_POINT:
	ar5 = ar1;
	ar2 = ar0 + 32;
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LANE
	FFT128_ROWX16_SCATTER_LAST
	ar0 = ar2;
	gr6 = gr6 + -1;
	gr6;
	if <>0 delayed goto LROWX16_SCATTER_POINT;
	ar1 = ar1 + 4;
	nul;
	gr5 = gr5 + -1;
	if <>0 delayed goto L28;
	ar6 = ar6 + 8192;
	nul;
	nul;
	ar1 = ar7 + -36;
	gr1 = [ar1];
	[ar7 + -12] = gr1;
	ar1 = ar7 + -35;
	gr1 = [ar1];
	[ar7 + -11] = gr1;
	ar0 = ar7 + -40;
	ar4 = [ar0];
	ar3 = _g_row_scratch;
	gr5 = 8;
L34:
	ar5 = ar3;
	ar6 = _g_col_batch;
	gr6 = 128;
	.align;
L30:
		ar0 = ar6 + 32;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar6++] = ar1,gr1;
		ar1,gr1 = [ar5++];
		[ar0++] = ar1,gr1;
		ar6 = ar0;
		ar5 = ar5 + 448;
	gr6 = gr6 - 1;
	if <>0 delayed goto L30;
	nul;
	nul;
		[_g_x16_col_out] = ar4;
		call _fft1d_128x16;
		ar4 = ar4 + 32;
		ar3 = ar3 + 64;
	gr5 = gr5 - 1;
	if <>0 delayed goto L34;
	nul;
	nul;
LRET_CONTEST:
	ar7 = ar7 + -64;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size	_contest_fft_2d, .-_contest_fft_2d
	.section	.bss.imu2,"aw",@nobits
	.align 3;
	.type	_g_tables_ready, @object
	.size	_g_tables_ready, 4
_g_tables_ready:
	.zero	4
	.align 3;
	.type	_g_x16_col_out, @object
	.size	_g_x16_col_out, 4
_g_x16_col_out:
	.zero	4
	.align 3;
	.type	_g_bitrev, @object
	.size	_g_bitrev, 512
_g_bitrev:
	.zero	512
	.align 3;
	.type	_g_column, @object
	.size	_g_column, 2048
_g_column:
	.zero	2048
	.align 3;
	.type	_g_column_f, @object
	.size	_g_column_f, 1024
_g_column_f:
	.zero	1024
	.align 3;
	.type	_g_fft_tmp, @object
	.size	_g_fft_tmp, 2048
_g_fft_tmp:
	.zero	2048
	.align 3;
	.type	_g_col_batch, @object
	.size	_g_col_batch, 32768
_g_col_batch:
	.zero	32768
	.section	.bss.emi,"aw",@nobits
	.align 3;
	.type	_g_row_scratch, @object
	.size	_g_row_scratch, 262144
_g_row_scratch:
	.zero	262144
	.section	.bss.imu4,"aw",@nobits
	.align 3;
	.type	_g_col_tmp, @object
	.size	_g_col_tmp, 32768
_g_col_tmp:
	.zero	32768
	.section	.bss.imu2,"aw",@nobits
	.align 3;
	.type	_g_roots, @object
	.size	_g_roots, 2048
_g_roots:
	.zero	2048
	.data
	.align 3;
	.type	_k_base_cos, @object
	.size	_k_base_cos, 512
_k_base_cos:
	.long	0
	.long	1072693248
	.long	3816385918
	.long	1072690721
	.long	2748392742
	.long	1072683149
	.long	2131272983
	.long	1072670549
	.long	3489094832
	.long	1072652951
	.long	4220661975
	.long	1072630398
	.long	1455828442
	.long	1072602945
	.long	83265253
	.long	1072570657
	.long	3476196678
	.long	1072533611
	.long	4097623923
	.long	1072491898
	.long	4051746225
	.long	1072445618
	.long	1159294526
	.long	1072394883
	.long	688824739
	.long	1072339814
	.long	1207142209
	.long	1072280544
	.long	1796544321
	.long	1072217216
	.long	577417135
	.long	1072149983
	.long	1719614413
	.long	1072079006
	.long	881650847
	.long	1072004457
	.long	621354454
	.long	1071926515
	.long	4242632758
	.long	1071845367
	.long	967731401
	.long	1071761211
	.long	2569207789
	.long	1071674247
	.long	939980350
	.long	1071524701
	.long	165764290
	.long	1071340816
	.long	2796464484
	.long	1071152610
	.long	1974149085
	.long	1070960538
	.long	785751813
	.long	1070765062
	.long	2065784606
	.long	1070537209
	.long	1013556749
	.long	1070135480
	.long	1854824762
	.long	1069729808
	.long	3156849720
	.long	1069094822
	.long	4044216357
	.long	1068048229
	.long	0
	.long	0
	.long	4044216339
	.long	-1079435419
	.long	3156849711
	.long	-1078388826
	.long	1854824758
	.long	-1077753840
	.long	1013556744
	.long	-1077348168
	.long	2065784602
	.long	-1076946439
	.long	785751810
	.long	-1076718586
	.long	1974149083
	.long	-1076523110
	.long	2796464482
	.long	-1076331038
	.long	165764284
	.long	-1076142832
	.long	939980348
	.long	-1075958947
	.long	2569207789
	.long	-1075809401
	.long	967731398
	.long	-1075722437
	.long	4242632757
	.long	-1075638281
	.long	621354453
	.long	-1075557133
	.long	881650848
	.long	-1075479191
	.long	1719614412
	.long	-1075404642
	.long	577417133
	.long	-1075333665
	.long	1796544321
	.long	-1075266432
	.long	1207142208
	.long	-1075203104
	.long	688824740
	.long	-1075143834
	.long	1159294525
	.long	-1075088765
	.long	4051746224
	.long	-1075038030
	.long	4097623923
	.long	-1074991750
	.long	3476196678
	.long	-1074950037
	.long	83265252
	.long	-1074912991
	.long	1455828442
	.long	-1074880703
	.long	4220661975
	.long	-1074853250
	.long	3489094832
	.long	-1074830697
	.long	2131272983
	.long	-1074813099
	.long	2748392741
	.long	-1074800499
	.long	3816385918
	.long	-1074792927
	.align 3;
	.type	_k_base_sin, @object
	.size	_k_base_sin, 512
_k_base_sin:
	.long	0
	.long	0
	.long	4044216341
	.long	1068048229
	.long	3156849708
	.long	1069094822
	.long	1854824762
	.long	1069729808
	.long	1013556746
	.long	1070135480
	.long	2065784602
	.long	1070537209
	.long	785751813
	.long	1070765062
	.long	1974149085
	.long	1070960538
	.long	2796464483
	.long	1071152610
	.long	165764288
	.long	1071340816
	.long	939980347
	.long	1071524701
	.long	2569207789
	.long	1071674247
	.long	967731400
	.long	1071761211
	.long	4242632757
	.long	1071845367
	.long	621354454
	.long	1071926515
	.long	881650847
	.long	1072004457
	.long	1719614412
	.long	1072079006
	.long	577417135
	.long	1072149983
	.long	1796544321
	.long	1072217216
	.long	1207142208
	.long	1072280544
	.long	688824739
	.long	1072339814
	.long	1159294526
	.long	1072394883
	.long	4051746224
	.long	1072445618
	.long	4097623923
	.long	1072491898
	.long	3476196678
	.long	1072533611
	.long	83265253
	.long	1072570657
	.long	1455828443
	.long	1072602945
	.long	4220661975
	.long	1072630398
	.long	3489094832
	.long	1072652951
	.long	2131272983
	.long	1072670549
	.long	2748392741
	.long	1072683149
	.long	3816385918
	.long	1072690721
	.long	0
	.long	1072693248
	.long	3816385918
	.long	1072690721
	.long	2748392742
	.long	1072683149
	.long	2131272983
	.long	1072670549
	.long	3489094832
	.long	1072652951
	.long	4220661975
	.long	1072630398
	.long	1455828443
	.long	1072602945
	.long	83265253
	.long	1072570657
	.long	3476196678
	.long	1072533611
	.long	4097623924
	.long	1072491898
	.long	4051746225
	.long	1072445618
	.long	1159294526
	.long	1072394883
	.long	688824741
	.long	1072339814
	.long	1207142209
	.long	1072280544
	.long	1796544322
	.long	1072217216
	.long	577417134
	.long	1072149983
	.long	1719614413
	.long	1072079006
	.long	881650849
	.long	1072004457
	.long	621354454
	.long	1071926515
	.long	4242632758
	.long	1071845367
	.long	967731400
	.long	1071761211
	.long	2569207790
	.long	1071674247
	.long	939980351
	.long	1071524701
	.long	165764287
	.long	1071340816
	.long	2796464485
	.long	1071152610
	.long	1974149090
	.long	1070960538
	.long	785751814
	.long	1070765062
	.long	2065784609
	.long	1070537209
	.long	1013556759
	.long	1070135480
	.long	1854824764
	.long	1069729808
	.long	3156849724
	.long	1069094822
	.long	4044216334
	.long	1068048229
	.ident	"GCC: (NMC SDK GCC ) 4.8.3"
