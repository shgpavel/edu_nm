	.file	"solution.c"
	.text
	.align 1;
	.macro FFT512_BUTTERFLY
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
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	.endm

	.macro FFT512X16_BUTTERFLY
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

	.macro FFT512X16_BUTTERFLY_PRELOAD
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

	.macro FFT512_COLX16_GATHER_LANE
	ar4,gr4 = [ar5++];
	[ar1++] = ar4,gr4;
	ar4,gr4 = [ar5++];
	[ar2++] = ar4,gr4;
	.endm

	.macro FFT512_ROWX16_GATHER_LANE
	ar4,gr4 = [ar5++];
	[ar1++] = ar4;
	[ar2++] = gr4;
	ar5 = ar5 + -1026;
	.endm

	.macro FFT512_ROWX16_SCATTER_LANE
	ar4,gr4 = [ar0++];
	[ar5++] = ar4,gr4;
	ar4,gr4 = [ar2++];
	[ar5++] = ar4,gr4;
	ar5 = ar5 + 2044;
	.endm

	.type	_fft1d_512, @function
_fft1d_512:
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
	ar2 = _g_bitrev;
	gr3 = _g_bitrev + 512;
	gr6 = _g_fft_tmp;
LFFT512_BITREV:
	gr0 = [ar2];
	gr0 = gr0 << 2;
	gr0 = gr6 + gr0;
	ar1 = gr0;
	ar4,gr4 = [ar0++];
	[ar1++] = ar4,gr4;
	ar4,gr4 = [ar0++];
	ar2 = ar2 + 1;
	gr0 = ar2;
	gr0 - gr3;
	if <>0 delayed goto LFFT512_BITREV;
	[ar1] = ar4,gr4;
	nul;
	gr6 = 0;
	vlen = gr6;
	ar0 = _g_fft_tmp;
	gr1 = 256;
	.align;
LFFT512_STAGE1:
	ar3 = ar0;
	ar4 = ar0 + 4;
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
	if <>0 delayed goto LFFT512_STAGE1;
	ar0 = ar0 + 8;
	nul;
	set fp_wait;
	set fp_branch;
	ar0 = _g_fft_tmp;
	gr1 = 128;
	.align;
LFFT512_STAGE2:
	ar3 = ar0;
	ar4 = ar0 + 8;
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
	ar3 = ar0 + 4;
	ar4 = ar0 + 12;
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
	if <>0 delayed goto LFFT512_STAGE2;
	ar0 = ar0 + 16;
	nul;
	set fp_wait;
	gr0 = 16;
	gr2 = 256;
	gr4 = 4;
	gr7 = 64;
	set fp_branch;
	.align;
LFFT512_STAGE:
	ar0 = _g_fft_tmp;
	gr3 = gr7;
	.align;
LFFT512_BLOCK_SETUP:
	ar1 = ar0 + gr0;
	ar3 = ar0;
	ar4 = ar1;
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
	ar0 = ar0 + 4;
	ar1 = ar1 + 4;
	ar2 = _g_roots;
	ar2 = ar2 + gr2;
	gr1 = gr4;
	gr1 = gr1 - 1;
	.align;
LFFT512_INNER:
	FFT512_BUTTERFLY
	ar0 = ar0 + 4;
	ar1 = ar1 + 4;
	ar2 = ar2 + gr2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512_INNER;
	nul;
	nul;
	gr3 = gr3 - 1;
	if =0 delayed goto LFFT512_STAGE_DONE;
	nul;
	nul;
	ar0 = ar1;
	delayed goto LFFT512_BLOCK_SETUP;
	nul;
	nul;
	.align;
LFFT512_STAGE_DONE:
	set fp_wait;
	gr0 = gr0 + gr0;
	gr2 = gr2 >> 1;
	gr4 = gr4 + gr4;
	gr7 = gr7 >> 1;
	gr7;
	if =0 delayed goto LFFT512_COPY_BACK_SETUP;
	nul;
	nul;
	gr6 = 1;
	gr7 - gr6;
	if =0 delayed goto LFFT512_STAGE9_DIRECT;
	nul;
	nul;
	set fp_branch;
	delayed goto LFFT512_STAGE;
	nul;
	nul;
LFFT512_STAGE9_DIRECT:
	set fp_branch;
	gr6 = [ar7 + -1];
	ar6 = gr6;
	ar0 = _g_fft_tmp;
	ar1 = _g_fft_tmp + 1024;
	ar4 = ar6 + 1024;
	gr6 = ar4;
	ar3 = ar0;
	ar4 = ar1;
	fpu 0 rep vlen vreg6 = [ar3++];
	fpu 0 rep vlen vreg7 = [ar3++];
	fpu 0 rep vlen vreg4 = [ar4++];
	fpu 0 rep vlen vreg5 = [ar4++];
	fpu 0 .double vreg0 = vreg6 + vreg4;
	fpu 0 .double vreg1 = vreg7 + vreg5;
	fpu 0 .double vreg2 = vreg6 - vreg4;
	fpu 0 .double vreg3 = vreg7 - vreg5;
	ar0 = ar0 + 4;
	ar1 = ar1 + 4;
	ar6 = ar6 + 4;
	ar4 = gr6;
	ar4 = ar4 + 4;
	gr6 = ar4;
	ar3 = ar6;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	ar2 = _g_roots + 4;
	gr1 = 255;
	.align;
LFFT512_STAGE9_INNER:
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
	ar0 = ar0 + 4;
	ar1 = ar1 + 4;
	ar2 = ar2 + 4;
	ar6 = ar6 + 4;
	ar4 = gr6;
	ar4 = ar4 + 4;
	gr6 = ar4;
	ar3 = ar6;
	fpu 0 rep vlen [--ar3] = vreg1;
	fpu 0 rep vlen [--ar3] = vreg0;
	fpu 0 rep vlen [--ar4] = vreg3;
	fpu 0 rep vlen [--ar4] = vreg2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512_STAGE9_INNER;
	nul;
	nul;
	set fp_wait;
	delayed goto LFFT512_EPILOGUE;
	nul;
	nul;
LFFT512_COPY_BACK_SETUP:
	set fp_wait;
	gr6 = [ar7 + -1];
	ar0 = gr6;
	ar1 = _g_fft_tmp;
	gr1 = 512;
LFFT512_COPY_BACK:
	ar5,gr5 = [ar1++];
	[ar0++] = ar5,gr5;
	ar5,gr5 = [ar1++];
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512_COPY_BACK;
	[ar0++] = ar5,gr5;
	nul;
LFFT512_EPILOGUE:
	ar7 = ar7 + -2;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size	_fft1d_512, .-_fft1d_512
	.align 1;
	.globl _fft1d_512x16
	.type	_fft1d_512x16, @function
_fft1d_512x16:
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
	gr1 = 128;
	.align;
LFFT512X16_BITREV_STAGE12:
	gr0 = [ar2++];
	gr0 = gr0 << 6;
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg0 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar3++];
	gr0 = gr0 << 6;
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg2 = [ar3++];
	fpu 0 rep vlen vreg3 = [ar3++];
	fpu 0 .double vreg4 = vreg0 + vreg2;
	fpu 0 .double vreg5 = vreg1 + vreg3;
	fpu 0 .double vreg6 = vreg0 - vreg2;
	fpu 0 .double vreg7 = vreg1 - vreg3;
	gr0 = gr0 << 6;
	gr0 = gr5 + gr0;
	ar3 = gr0;
	gr0 = [ar2++];
	fpu 0 rep vlen vreg0 = [ar3++];
	fpu 0 rep vlen vreg1 = [ar3++];
	gr0 = gr0 << 6;
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
	if <>0 delayed goto LFFT512X16_BITREV_STAGE12;
	ar1 = ar1 + 192;
	nul;
	set fp_wait;
	set fp_branch;
	ar0 = _g_col_tmp;
	gr1 = 64;
	.align;
LFFT512X16_STAGE3_R0:
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
	if <>0 delayed goto LFFT512X16_STAGE3_R0;
	ar0 = ar0 + 512;
	nul;
	set fp_wait;
	ar2 = _g_roots + 256;
	gr1 = 64;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 64;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 320;
	set fp_branch;
	.align;
LFFT512X16_STAGE3_R64:
	FFT512X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 512;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_STAGE3_R64;
	ar0 = ar0 + 512;
	nul;
	set fp_wait;
	ar0 = _g_col_tmp + 128;
	gr1 = 64;
	set fp_branch;
	.align;
LFFT512X16_STAGE3_R128:
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
	if <>0 delayed goto LFFT512X16_STAGE3_R128;
	ar0 = ar0 + 512;
	nul;
	set fp_wait;
	ar2 = _g_roots + 768;
	gr1 = 64;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 192;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 448;
	set fp_branch;
	.align;
LFFT512X16_STAGE3_R192:
	FFT512X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 512;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_STAGE3_R192;
	ar0 = ar0 + 512;
	nul;
	set fp_wait;
	set fp_branch;
	ar0 = _g_col_tmp;
	gr1 = 32;
	.align;
LFFT512X16_STAGE4_R0:
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
	if <>0 delayed goto LFFT512X16_STAGE4_R0;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 128;
	gr1 = 32;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 64;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 576;
	set fp_branch;
	.align;
LFFT512X16_STAGE4_R32:
	FFT512X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_STAGE4_R32;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 256;
	gr1 = 32;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 128;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 640;
	set fp_branch;
	.align;
LFFT512X16_STAGE4_R64:
	FFT512X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_STAGE4_R64;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 384;
	gr1 = 32;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 192;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 704;
	set fp_branch;
	.align;
LFFT512X16_STAGE4_R96:
	FFT512X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_STAGE4_R96;
	ar0 = ar0 + 1024;
	nul;
	ar0 = _g_col_tmp + 256;
	gr1 = 32;
	set fp_branch;
	.align;
LFFT512X16_STAGE4_R128:
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
	if <>0 delayed goto LFFT512X16_STAGE4_R128;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 640;
	gr1 = 32;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 320;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 832;
	set fp_branch;
	.align;
LFFT512X16_STAGE4_R160:
	FFT512X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_STAGE4_R160;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 768;
	gr1 = 32;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 384;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 896;
	set fp_branch;
	.align;
LFFT512X16_STAGE4_R192:
	FFT512X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_STAGE4_R192;
	ar0 = ar0 + 1024;
	nul;
	ar2 = _g_roots + 896;
	gr1 = 32;
	fpu 0 rep vlen vreg4 = [ar2];
	ar2 = ar2 + 2;
	ar0 = _g_col_tmp + 448;
	fpu 0 rep vlen vreg5 = [ar2];
	ar1 = _g_col_tmp + 960;
	set fp_branch;
	.align;
LFFT512X16_STAGE4_R224:
	FFT512X16_BUTTERFLY_PRELOAD
	ar1 = ar1 + 1024;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_STAGE4_R224;
	ar0 = ar0 + 1024;
	nul;
	set fp_wait;
	gr0 = 1024;
	gr2 = 64;
	gr4 = 16;
	gr7 = 16;
	set fp_branch;
	.align;
LFFT512X16_STAGE:
	ar0 = _g_col_tmp;
	gr3 = gr7;
	.align;
LFFT512X16_BLOCK_SETUP:
	ar1 = ar0 + gr0;
	ar3 = ar0;
	ar4 = ar1;
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
	ar0 = ar0 + 64;
	ar1 = ar1 + 64;
	ar2 = _g_roots;
	ar2 = ar2 + gr2;
	gr1 = gr4;
	gr1 = gr1 - 1;
	.align;
LFFT512X16_INNER:
	FFT512X16_BUTTERFLY
	ar0 = ar0 + 64;
	ar1 = ar1 + 64;
	ar2 = ar2 + gr2;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_INNER;
	nul;
	nul;
	gr3 = gr3 - 1;
	if =0 delayed goto LFFT512X16_STAGE_DONE;
	nul;
	nul;
	ar0 = ar1;
	delayed goto LFFT512X16_BLOCK_SETUP;
	nul;
	nul;
	.align;
LFFT512X16_STAGE_DONE:
	set fp_wait;
	gr0 = gr0 + gr0;
	gr2 = gr2 >> 1;
	gr4 = gr4 + gr4;
	gr7 = gr7 >> 1;
	gr7;
	if =0 delayed goto LFFT512X16_RET;
	gr6 = gr7 + -1;
	set fp_branch;
	if <>0 delayed goto LFFT512X16_STAGE;
	ar5 = [_g_x16_col_out];
	ar0 = _g_col_tmp;
	gr6 = ar5;
	gr6;
	if =0 delayed goto LFFT512X16_STAGE;
	ar1 = ar0 + gr0;
	ar2 = _g_roots;
	.align;
LFFT512X16_COL_FINAL:
	gr6 = 15;
	vlen = gr6;
	set fp_branch;
	sir = 0;
	fp_pack_exp = sir;
	ar6 = ar5 + 262144;
	gr1 = 256;
	.align;
LFFT512X16_COL_FINAL_INNER:
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
	fpu 0 .packer = (vreg0,vreg1) with .fixed_32 <= .double;
	ar3 = _g_fft_tmp;
	fpu rep 16 [ar3++] = .packer;
	ar3 = _g_fft_tmp;
	fpu rep 16 .packer = [ar3++] with .float <= .fixed_32;
	fpu rep 16 [ar5++] = .packer;
	ar5 = ar5 + 992;
	fpu 0 .packer = (vreg2,vreg3) with .fixed_32 <= .double;
	ar3 = _g_fft_tmp;
	fpu rep 16 [ar3++] = .packer;
	ar3 = _g_fft_tmp;
	fpu rep 16 .packer = [ar3++] with .float <= .fixed_32;
	fpu rep 16 [ar6++] = .packer;
	gr1 = gr1 - 1;
	if <>0 delayed goto LFFT512X16_COL_FINAL_INNER;
	ar6 = ar6 + 992;
	nul;
	[_g_x16_col_out] = gr1;
	set fp_wait;
	nul;
LFFT512X16_RET:
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size	_fft1d_512x16, .-_fft1d_512x16
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
	ar7 = ar7 + 160;
	ar6 = ar7 + -160;
	ar0 = ar6 + -16;
	gr7 = [ar0 + -1];
	[ar6 + 80] = gr7;
	ar2 = [_g_tables_ready];
	gr6 = ar2;
	gr6;
	if <>0 goto L15;
	ar0 = _g_roots;
	ar1 = _k_quarter_sin + 256;
	ar3 = _k_quarter_sin;
	[ar6 + 84] = gr4;
	[ar6 + 85] = gr5;
	gr4 = ar1;
	ar2 = ar6 + 140;
	gr7 = [ar2];
	[ar6 + 100] = gr7;
	ar2 = ar6 + 141;
	gr7 = [ar2];
	[ar6 + 101] = gr7;
	gr5 = gr6;
	ar5 = ar0;
	gr6 = ar3;
	[ar6 + 92] = gr0;
	[ar6 + 93] = gr1;
L20:
	gr7 = 128;
	gr7 - gr5;
	if u< goto L16;
	ar0 = gr4;
	ar2 = [ar0];
	ar3 = [ar0 + 1];
	ar1 = gr6;
	ar0 = [ar1];
	[ar6 + 148] = ar0;
	ar0 = [ar1 + 1];
	[ar6 + 149] = ar0;
	[ar5] = ar2;
	[ar5 + 1] = ar3;
	ar1 = ar5 + 2;
	ar4 = [ar6 + 148];
	[ar1] = ar4;
	ar0 = [ar6 + 149];
	[ar5 + 3] = ar0;
	ar1 = ar5 + 1024;
	ar4 = ar3 + -2147483648;
	[ar1] = ar2;
	[ar5 + 1025] = ar4;
	ar4 = ar5 + 1026;
	[ar6 + 108] = ar4;
	ar0 = ar6 + 148;
	gr2 = [ar0];
	ar4 = [ar6 + 149];
	ar1 = ar4 + -2147483648;
	gr3 = ar1;
	ar1 = [ar6 + 108];
	[ar1] = gr2;
	[ar5 + 1027] = gr3;
	ar0 = gr5;
	ar0 = ar0 + 1;
	gr5 = ar0;
	ar5 = ar5 + 4;
	ar1 = gr4;
	ar1 = ar1 + -2;
	gr4 = ar1;
	ar4 = gr6;
	ar4 = ar4 + 2;
	gr6 = ar4;
	gr7 = 256;
	gr5 - gr7;
	if <>0 goto L20;
L47:
	ar0 = ar6 + 84;
	gr4 = [ar0];
	ar1 = ar6 + 85;
	gr5 = [ar1];
	ar3 = ar6 + 100;
	gr1 = [ar3];
	[ar6 + 140] = gr1;
	ar1 = ar6 + 101;
	gr1 = [ar1];
	[ar6 + 141] = gr1;
	ar0 = ar6 + 92;
	gr0 = [ar0];
	ar1 = ar6 + 93;
	gr1 = [ar1];
	ar1 = 0;
	ar2 = _g_bitrev;
	ar3 = 1;
	[ar6 + 72] = gr0;
	[ar6 + 73] = gr1;
L19:
	ar0 = 9;
	gr0 = 0;
	gr6 = gr0;
	gr1 = ar1;
L23:
	gr2 = gr6 + gr6;
	gr6 = ar3;
	gr0 = gr1 and gr6;
	gr6 = gr0 or gr2;
	gr1 = gr1 >> 1;
	ar0 = ar0 + -1;
	gr7 = ar0;
	gr7;
	if <>0 goto L23;
	gr0 = gr6;
	gr6 = ar2;
	gr7 = ar1;
	gr1 = gr6 + gr7;
	ar0 = gr1;
	[ar0] = gr0;
	ar1 = ar1 + 1;
	gr6 = 512;
	gr7 = ar1;
	gr7 - gr6;
	if <>0 goto L19;
	gr6 = 1;
	[_g_tables_ready] = gr6;
L15:
	ar3 = ar6 + 80;
	gr3 = [ar3];
	[ar6 + 72] = gr3;
	gr6 = _g_row_scratch;
	ar2 = ar7 + -1;
	[ar6 + 84] = ar2;
	[ar6 + 100] = gr4;
	[ar6 + 101] = gr5;
	gr5 = gr6;
	ar1 = ar6 + 140;
	gr3 = [ar1];
	[ar6 + 92] = gr3;
	ar1 = ar6 + 141;
	gr3 = [ar1];
	[ar6 + 93] = gr3;
	set fp_branch;
	ar3 = [ar6 + 72];
	gr5 = _g_row_scratch;
	gr7 = 32;
L27X16:
	[ar6 + 72] = ar3;
	[ar6 + 76] = gr5;
	ar0 = _g_col_batch;
	gr6 = 512;
	.align;
L25X16:
	ar5 = ar3;
	ar5 = ar5 + 15360;
	ar1 = _g_fft_tmp;
	ar2 = _g_fft_tmp + 16;
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	FFT512_ROWX16_GATHER_LANE
	ar1 = _g_fft_tmp;
	ar2 = ar0;
	fpu rep 8 .packer = [ar1++] with .double <= .float;
	fpu rep 16 [ar2++] = .packer;
	ar1 = _g_fft_tmp + 16;
	fpu rep 8 .packer = [ar1++] with .double <= .float;
	fpu rep 16 [ar2++] = .packer;
	ar0 = ar2;
	gr6 = gr6 - 1;
	if <>0 delayed goto L25X16;
	ar3 = ar3 + 2;
	nul;
	set fp_wait;
	gr6 = 0;
	[_g_x16_col_out] = gr6;
	[ar6 + 88] = gr7;
	call _fft1d_512x16;
	ar1 = ar6 + 88;
	gr7 = [ar1];
	ar3 = [ar6 + 76];
	ar0 = _g_col_tmp;
	gr6 = 512;
	.align;
LROWX16_SCATTER_POINT:
	ar5 = ar3;
	ar2 = ar0 + 32;
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	FFT512_ROWX16_SCATTER_LANE
	ar0 = ar2;
	gr6 = gr6 - 1;
	if <>0 delayed goto LROWX16_SCATTER_POINT;
	ar3 = ar3 + 4;
	nul;
	ar3 = [ar6 + 72];
	ar3 = ar3 + 16384;
	gr5 = [ar6 + 76];
	ar0 = gr5;
	ar0 = ar0 + 32768;
	gr5 = ar0;
	gr7 = gr7 - 1;
	if <>0 delayed goto L27X16;
	nul;
	nul;
	ar3 = ar6 + 92;
	gr1 = [ar3];
	[ar6 + 140] = gr1;
	ar1 = ar6 + 93;
	gr1 = [ar1];
	[ar6 + 141] = gr1;
	ar3 = _g_row_scratch;
	gr1 = _g_column;
	[ar6 + 72] = gr1;
	ar3 = _g_row_scratch;
	ar4 = [ar6 + 80];
	gr5 = 32;
L39X16:
	[ar6 + 72] = ar3;
	[ar6 + 76] = ar4;
	ar0 = _g_col_batch;
	gr6 = 512;
	.align;
L29X16:
	ar5 = ar3;
	ar1 = ar0;
	ar2 = ar0 + 32;
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	FFT512_COLX16_GATHER_LANE
	ar0 = ar0 + 64;
	gr6 = gr6 - 1;
	if <>0 delayed goto L29X16;
	ar3 = ar3 + 2048;
	nul;
	set fp_wait;
	ar4 = [ar6 + 76];
	[_g_x16_col_out] = ar4;
	call _fft1d_512x16;
	ar4 = ar4 + 32;
	ar3 = [ar6 + 72];
	ar3 = ar3 + 64;
	gr5 = gr5 - 1;
	if <>0 delayed goto L39X16;
	nul;
	nul;
	goto LRET_CONTEST;
L39:
	ar0 = ar3;
	ar2 = ar6 + 80;
	ar1 = [ar2];
	ar2 = _g_column;
	ar4 = ar3 + 1048576;
	gr7 = ar4;
L29:
	ar5,gr5 = [ar0++];
	[ar2++] = ar5,gr5;
	ar5,gr5 = [ar0++];
	[ar2++] = ar5,gr5;
	ar0 = ar0 + 2044;
	gr0 = ar0;
	gr0 - gr7;
	if <>0 delayed goto L29;
	nul;
	nul;
	gr1 = [ar6 + 72];
	ar2 = [ar6 + 84];
	[ar2] = gr1;
	call _fft1d_512;
	gr1 = _g_column;
	[ar6 + 72] = gr1;
	ar0 = gr1;
	[ar6 + 92] = ar3;
L37:
	set fp_branch;
	sir = 0;
	fp_pack_exp = sir;
	ar0 = _g_column;
	ar2 = _g_fft_tmp;
	gr2 = 32;
LROUND_TO_FIXED:
	fpu 0 rep 32 vreg0 = [ar0++];
	fpu 0 .packer = vreg0 with .fixed_32 <= .double;
	fpu rep 16 [ar2++] = .packer;
	gr2 = gr2 - 1;
	if <>0 delayed goto LROUND_TO_FIXED;
	nul;
	nul;
	ar0 = _g_fft_tmp;
	gr7 = _g_fft_tmp + 1024;
	ar2 = gr7;
	gr2 = 32;
LFIXED_TO_FLOAT:
	fpu rep 16 .packer = [ar0++] with .float <= .fixed_32;
	fpu rep 16 [ar2++] = .packer;
	gr2 = gr2 - 1;
	if <>0 delayed goto LFIXED_TO_FLOAT;
	nul;
	nul;
	set fp_wait;
	gr7 = _g_fft_tmp + 1024;
	ar0 = gr7;
	ar1 = [ar6 + 80];
	gr2 = 512;
LSCATTER_COLUMN:
	ar5,gr5 = [ar0++];
	[ar1++] = ar5,gr5;
	gr2 = gr2 - 1;
	if <>0 delayed goto LSCATTER_COLUMN;
	ar1 = ar1 + 1022;
	nul;
L46:
	ar0 = ar6 + 92;
	ar3 = [ar0];
	ar1 = [ar6 + 80];
	ar4 = ar1 + 2;
	[ar6 + 80] = ar4;
	ar3 = ar3 + 4;
	gr1 = ar3;
	gr2 = _g_row_scratch + 2048;
	gr1 - gr2;
	if <>0 goto L39;
LRET_CONTEST:
	ar7 = ar7 + -160;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L16:
	ar1 = gr6;
	ar0 = [ar1 + -256];
	[ar6 + 72] = ar0;
	ar0 = [ar1 + -255];
	[ar6 + 108] = ar0;
	ar3 = ar0 + -2147483648;
	[ar6 + 73] = ar3;
	ar0 = ar6 + 72;
	ar2 = [ar0];
	ar4 = gr4;
	ar1 = [ar4 + 256];
	[ar6 + 148] = ar1;
	ar1 = [ar4 + 257];
	[ar6 + 149] = ar1;
	[ar5] = ar2;
	[ar5 + 1] = ar3;
	ar1 = ar5 + 2;
	ar4 = [ar6 + 148];
	[ar1] = ar4;
	ar0 = [ar6 + 149];
	[ar5 + 3] = ar0;
	ar1 = ar5 + 1024;
	ar4 = ar3 + -2147483648;
	[ar1] = ar2;
	[ar5 + 1025] = ar4;
	ar4 = ar5 + 1026;
	[ar6 + 108] = ar4;
	ar0 = ar6 + 148;
	gr2 = [ar0];
	ar4 = [ar6 + 149];
	ar1 = ar4 + -2147483648;
	gr3 = ar1;
	ar1 = [ar6 + 108];
	[ar1] = gr2;
	[ar5 + 1027] = gr3;
	ar0 = gr5;
	ar0 = ar0 + 1;
	gr5 = ar0;
	ar5 = ar5 + 4;
	ar1 = gr4;
	ar1 = ar1 + -2;
	gr4 = ar1;
	ar4 = gr6;
	ar4 = ar4 + 2;
	gr6 = ar4;
	gr7 = 256;
	gr5 - gr7;
	if <>0 goto L20;
	goto L47;
	.size	_contest_fft_2d, .-_contest_fft_2d
	.section	.bss.emi,"aw",@nobits
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
	.size	_g_bitrev, 2048
_g_bitrev:
	.zero	2048
	.align 3;
	.type	_g_row_scratch, @object
	.size	_g_row_scratch, 4194304
_g_row_scratch:
	.zero	4194304
	.align 3;
	.type	_g_column, @object
	.size	_g_column, 8192
_g_column:
	.zero	8192
	.align 3;
	.type	_g_fft_tmp, @object
	.size	_g_fft_tmp, 8192
_g_fft_tmp:
	.zero	8192
	.align 3;
	.type	_g_col_batch, @object
	.size	_g_col_batch, 131072
_g_col_batch:
	.zero	131072
	.align 3;
	.type	_g_col_tmp, @object
	.size	_g_col_tmp, 131072
_g_col_tmp:
	.zero	131072
	.align 3;
	.type	_g_roots, @object
	.size	_g_roots, 8192
_g_roots:
	.zero	8192
	.data
	.align 3;
	.type	_k_quarter_sin, @object
	.size	_k_quarter_sin, 1032
_k_quarter_sin:
	.long	0
	.long	0
	.long	4242458500
	.long	1065951697
	.long	4154680958
	.long	1067000149
	.long	1972655565
	.long	1067636837
	.long	4044216340
	.long	1068048229
	.long	2046304480
	.long	1068459374
	.long	2462980598
	.long	1068684576
	.long	177906713
	.long	1068889808
	.long	3156849708
	.long	1069094822
	.long	3349070549
	.long	1069299589
	.long	1453945614
	.long	1069504078
	.long	3470456182
	.long	1069627888
	.long	1854824762
	.long	1069729808
	.long	2972270170
	.long	1069831542
	.long	1149976518
	.long	1069933076
	.long	3729833777
	.long	1070034393
	.long	1013556746
	.long	1070135480
	.long	632292433
	.long	1070236320
	.long	1786662755
	.long	1070336898
	.long	3846521617
	.long	1070437199
	.long	2065784602
	.long	1070537209
	.long	238550446
	.long	1070616504
	.long	1655830135
	.long	1070666194
	.long	1062001469
	.long	1070715716
	.long	785751813
	.long	1070765062
	.long	3269558610
	.long	1070814224
	.long	2484561410
	.long	1070863196
	.long	1115219818
	.long	1070911970
	.long	1974149085
	.long	1070960538
	.long	3711903686
	.long	1071008893
	.long	821708190
	.long	1071057029
	.long	824036583
	.long	1071104937
	.long	2796464483
	.long	1071152610
	.long	1673302999
	.long	1071200042
	.long	1135144330
	.long	1071247225
	.long	728580041
	.long	1071294152
	.long	165764288
	.long	1071340816
	.long	3623952103
	.long	1071387209
	.long	2570175581
	.long	1071433326
	.long	1535566729
	.long	1071479159
	.long	939980347
	.long	1071524701
	.long	1391425749
	.long	1071569945
	.long	3690502842
	.long	1071614884
	.long	2269920951
	.long	1071652092
	.long	2569207789
	.long	1071674247
	.long	1297975695
	.long	1071696240
	.long	1411074850
	.long	1071718067
	.long	1675549513
	.long	1071739725
	.long	967731400
	.long	1071761211
	.long	2570316632
	.long	1071782521
	.long	1289556703
	.long	1071803653
	.long	637203434
	.long	1071824603
	.long	4242632757
	.long	1071845367
	.long	2969983475
	.long	1071865944
	.long	689983673
	.long	1071886330
	.long	1397053140
	.long	1071906521
	.long	621354454
	.long	1071926515
	.long	2610628920
	.long	1071946308
	.long	3152274921
	.long	1071965898
	.long	2460177814
	.long	1071985282
	.long	881650847
	.long	1072004457
	.long	3194290572
	.long	1072023419
	.long	1427975390
	.long	1072042167
	.long	636516414
	.long	1072060697
	.long	1719614412
	.long	1072079006
	.long	1424664750
	.long	1072097092
	.long	938475414
	.long	1072114952
	.long	1594061409
	.long	1072132583
	.long	577417135
	.long	1072149983
	.long	3814135664
	.long	1072167148
	.long	201299821
	.long	1072184078
	.long	263859903
	.long	1072200768
	.long	1796544321
	.long	1072217216
	.long	2750387212
	.long	1072233420
	.long	1234330619
	.long	1072249378
	.long	4106737473
	.long	1072265086
	.long	1207142208
	.long	1072280544
	.long	4012486548
	.long	1072295747
	.long	2688952194
	.long	1072310695
	.long	453180130
	.long	1072325385
	.long	688824739
	.long	1072339814
	.long	2653017360
	.long	1072353980
	.long	1477771775
	.long	1072367882
	.long	1056266002
	.long	1072381517
	.long	1159294526
	.long	1072394883
	.long	1731563730
	.long	1072407978
	.long	2892993833
	.long	1072420800
	.long	645027145
	.long	1072433348
	.long	4051746224
	.long	1072445618
	.long	881357723
	.long	1072457611
	.long	557060597
	.long	1072469323
	.long	4093443164
	.long	1072480752
	.long	4097623923
	.long	1072491898
	.long	1950234075
	.long	1072502759
	.long	3511535930
	.long	1072513332
	.long	1942611595
	.long	1072523617
	.long	3476196678
	.long	1072533611
	.long	1647878298
	.long	1072543314
	.long	3066872380
	.long	1072552723
	.long	1647164971
	.long	1072561838
	.long	83265253
	.long	1072570657
	.long	1261158782
	.long	1072579178
	.long	3964199186
	.long	1072587400
	.long	2873938188
	.long	1072595323
	.long	1455828443
	.long	1072602945
	.long	3370060660
	.long	1072610264
	.long	3882371427
	.long	1072617280
	.long	2749658050
	.long	1072623992
	.long	4220661975
	.long	1072630398
	.long	4151720665
	.long	1072636498
	.long	2892293513
	.long	1072642291
	.long	990588717
	.long	1072647776
	.long	3489094832
	.long	1072652951
	.long	2745245815
	.long	1072657817
	.long	4201728937
	.long	1072662372
	.long	617154971
	.long	1072666617
	.long	2131272983
	.long	1072670549
	.long	1905645377
	.long	1072674169
	.long	1893834707
	.long	1072677476
	.long	4251821691
	.long	1072680469
	.long	2748392741
	.long	1072683149
	.long	4240268041
	.long	1072685514
	.long	2902558778
	.long	1072687565
	.long	1703834029
	.long	1072689301
	.long	3816385918
	.long	1072690721
	.long	4026464093
	.long	1072691826
	.long	1619315981
	.long	1072692616
	.long	379294427
	.long	1072693090
	.long	0
	.long	1072693248
	.ident	"GCC: (NMC SDK GCC ) 4.8.3"
