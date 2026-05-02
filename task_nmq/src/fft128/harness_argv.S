	.file	"harness_argv.c"
	.text
	.align 1;
	.type	_fft1d_reference, @function
_fft1d_reference:
	push ar6,gr6;
	push ar5,gr5;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	ar7 = ar7 + 68;
	ar6 = ar7 + -68;
	ar0 = ar6 + -16;
	gr1 = [ar0 + -1];
	[ar6 + 48] = gr1;
	gr0 = 0;
	gr3 = _vec_buf0;
	gr2 = 256;
	gr4 = gr1;
L4:
	gr1 = gr3 + gr0;
	ar0 = gr1;
	gr5 = gr4 + gr0;
	ar1 = gr5;
	gr1 = [ar1];
	[ar0] = gr1;
	ar1 = [ar1 + 1];
	[ar0 + 1] = ar1;
	ar0 = gr0;
	ar0 = ar0 + 2;
	gr0 = ar0;
	gr0 - gr2;
	if <>0 goto L4;
	gr1 = 7;
	[ar6 + 40] = gr1;
	gr1 = 64;
	[ar6 + 8] = gr1;
	gr1 = 1;
	[ar6 + 32] = gr1;
	gr1 = _vec_buf1;
	[ar6 + 36] = gr1;
	gr1 = _vec_buf0;
	[ar6 + 44] = gr1;
L3:
	gr0 = [ar6 + 32];
	gr0;
	if v<= goto L9;
	ar1 = [ar6 + 44];
	ar0 = ar1 + 1;
	gr1 = ar0;
	ar0 = ar6 + 36;
	gr0 = [ar0];
	gr3 = [ar6 + 8];
	gr3 = gr3 + gr3;
	[ar6 + 16] = gr3;
	gr3 = 0;
	[ar6 + 12] = gr3;
	[ar6 + 20] = gr3;
L10:
	gr3 = _twiddle_cos;
	gr4 = [ar6 + 12];
	gr2 = gr3 + gr4;
	ar0 = gr2;
	ar0 = [ar0];
	[ar6] = ar0;
	gr3 = _twiddle_sin;
	gr2 = gr3 + gr4;
	ar1 = gr2;
	ar0 = [ar1];
	ar4 = ar0 + -2147483648;
	[ar6 + 4] = ar4;
	gr3 = [ar6 + 16];
	gr3 = gr1 + gr3;
	[ar6 + 28] = gr3;
	gr3 = [ar6 + 16];
	gr3 = gr0 + gr3;
	[ar6 + 24] = gr3;
	ar3 = gr3;
	ar2 = gr0;
	ar4 = ar6 + 28;
	ar1 = [ar4];
	ar0 = gr1;
	gr5 = 0;
	[ar6 + 64] = gr5;
	[ar6 + 60] = ar0;
L6:
	ar0 = [ar6 + 60];
	gr1 = [ar0 + -1];
	gr0 = [ar0];
	gr3 = [ar1 + -1];
	gr4 = [ar6];
	delayed call FMul;
	push gr4;
	push gr3;
	gr4 = gr7;
	[ar6 + 52] = gr4;
	gr4 = [ar6 + 4];
	gr2 = [ar1];
	delayed call FMul;
	push gr4;
	push gr2;
	gr4 = gr7;
	[ar6 + 56] = gr4;
	gr2 = [ar6 + 52];
	gr7 = 80000000h;
	gr7 = gr7 xor gr4;
	delayed call FAdd;
	push gr2;
	push gr7;
	gr4 = gr7;
	gr5 = [ar6 + 4];
	delayed call FMul;
	push gr5;
	push gr3;
	gr5 = gr7;
	[ar6 + 52] = gr5;
	gr3 = [ar6];
	gr2 = [ar1];
	delayed call FMul;
	push gr3;
	push gr2;
	gr3 = gr7;
	[ar6 + 56] = gr3;
	delayed call FAdd;
	push gr5;
	push gr3;
	gr2 = gr7;
	delayed call FAdd;
	push gr4;
	push gr1;
	gr3 = gr7;
	[ar2] = gr3;
	delayed call FAdd;
	push gr2;
	push gr0;
	gr3 = gr7;
	[ar2 + 1] = gr3;
	gr7 = 80000000h;
	gr7 = gr7 xor gr4;
	delayed call FAdd;
	push gr1;
	push gr7;
	gr4 = gr7;
	[ar3] = gr4;
	gr7 = 80000000h;
	gr7 = gr7 xor gr2;
	delayed call FAdd;
	push gr0;
	push gr7;
	gr0 = gr7;
	[ar3 + 1] = gr0;
	ar0 = [ar6 + 64];
	ar5 = ar0 + 1;
	[ar6 + 64] = ar5;
	ar0 = [ar6 + 60];
	ar5 = ar0 + 2;
	[ar6 + 60] = ar5;
	ar1 = ar1 + 2;
	ar2 = ar2 + 2;
	ar3 = ar3 + 2;
	gr0 = [ar6 + 64];
	gr1 = [ar6 + 8];
	gr0 - gr1;
	if v< goto L6;
	ar0 = [ar6 + 20];
	ar2 = ar0 + 1;
	[ar6 + 20] = ar2;
	gr0 = [ar6 + 12];
	gr1 = gr0 + gr1;
	[ar6 + 12] = gr1;
	gr0 = [ar6 + 28];
	gr2 = [ar6 + 16];
	gr1 = gr0 + gr2;
	gr3 = [ar6 + 24];
	gr0 = gr3 + gr2;
	gr2 = ar2;
	gr3 = [ar6 + 32];
	gr2 - gr3;
	if <>0 goto L10;
L9:
	gr1 = [ar6 + 32];
	gr1 = gr1 + gr1;
	[ar6 + 32] = gr1;
	gr1 = [ar6 + 8];
	gr1 = gr1 A>> 1;
	[ar6 + 8] = gr1;
	ar0 = [ar6 + 40];
	ar2 = ar0 + -1;
	[ar6 + 40] = ar2;
	gr0 = ar2;
	gr0;
	if =0 goto L8;
	ar0 = ar6 + 44;
	gr0 = [ar0];
	ar1 = ar6 + 36;
	gr2 = [ar1];
	[ar6 + 44] = gr2;
	[ar6 + 36] = gr0;
	goto L3;
L8:
	ar0 = ar6 + 48;
	ar2 = [ar0];
	ar3 = ar6 + 40;
	ar1 = [ar3];
	gr3 = 1;
	ar4 = ar6 + 36;
	gr4 = [ar4];
L14:
	gr1 = ar1;
	ar0 = 7;
	gr0 = 0;
L12:
	gr0 = gr0 + gr0;
	gr2 = gr1 and gr3;
	gr0 = gr2 or gr0;
	gr1 = gr1 A>> 1;
	ar0 = ar0 + -1;
	gr2 = ar0;
	gr2;
	if <>0 goto L12;
	gr0 = gr0 + gr0;
	gr0 = gr4 + gr0;
	ar0 = gr0;
	gr0 = [ar0];
	[ar2] = gr0;
	ar0 = [ar0 + 1];
	[ar2 + 1] = ar0;
	ar1 = ar1 + 1;
	ar2 = ar2 + 2;
	gr0 = 128;
	gr5 = ar1;
	gr5 - gr0;
	if <>0 goto L14;
	ar7 = ar7 + -68;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size	_fft1d_reference, .-_fft1d_reference
	.align 1;
	.type	_fill_case.constprop.2, @function
_fill_case.constprop.2:
	push ar6,gr6;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	ar6 = ar7;
	ar0 = ar6 + -12;
	gr2 = [ar0 + -1];
	ar0 = _template_buf;
	ar1 = ar0 + 32768;
	gr1 = ar1;
	gr0 = 000000000h;
L19:
	[ar0] = gr0;
	[ar0 + 1] = gr0;
	ar0 = ar0 + 2;
	gr3 = ar0;
	gr3 - gr1;
	if <>0 goto L19;
	gr2;
	if =0 goto L50;
	gr0 = 1;
	gr2 - gr0;
	if =0 goto L51;
	gr3 = 2;
	gr2 - gr3;
	if =0 goto L52;
	gr3 = 3;
	gr2 - gr3;
	if =0 goto L53;
	gr3 = 4;
	gr2 - gr3;
	if =0 goto L54;
	gr3 = 6;
	gr2 - gr3;
	if =0 goto L55;
	gr3 = 7;
	gr2 - gr3;
	if =0 goto L56;
	gr3 = 5;
	gr2 - gr3;
	if =0 goto L57;
	gr7 = gr0;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar6,gr6;
	return;
L52:
	ar1 = _template_buf;
	gr0 = 0;
L27:
	ar0 = ar1;
	gr6 = 0;
	gr1 = 1;
L31:
	gr7 = gr6 + gr0;
	gr7 = gr7 and gr1;
	gr7;
	if =0 goto L47;
	gr7 = 0bf800000h;
	[ar0] = gr7;
	ar2 = gr6;
	ar2 = ar2 + 1;
	gr6 = ar2;
	ar0 = ar0 + 2;
	gr7 = 128;
	gr6 - gr7;
	if <>0 goto L31;
L58:
	ar3 = gr0;
	ar3 = ar3 + 1;
	gr0 = ar3;
	ar1 = ar1 + 256;
	gr0 - gr6;
	if <>0 goto L27;
L46:
	gr7 = 0;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar6,gr6;
	return;
L47:
	gr7 = 03f800000h;
	[ar0] = gr7;
	ar2 = gr6;
	ar2 = ar2 + 1;
	gr6 = ar2;
	ar0 = ar0 + 2;
	gr7 = 128;
	gr6 - gr7;
	if <>0 goto L31;
	goto L58;
L50:
	gr0 = 03f800000h;
	[_template_buf] = gr0;
	gr7 = gr2;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar6,gr6;
	return;
L51:
	ar1 = _template_buf;
	gr6 = 03f800000h;
	gr0 = 128;
	ar0 = ar1;
	ar5 = 0;
L25:
	[ar0] = gr6;
	ar5 = ar5 + 1;
	ar0 = ar0 + 2;
	gr2 = ar5;
	gr2 - gr0;
	if <>0 goto L25;
	ar1 = ar1 + 256;
	gr3 = ar1;
	gr3 - gr1;
	if =0 goto L46;
	ar0 = ar1;
	ar5 = 0;
	goto L25;
L53:
	gr0 = 03f800000h;
	[_template_buf + 8192] = gr0;
	gr7 = 0;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar6,gr6;
	return;
L57:
	ar2 = _template_buf;
	gr2 = 15;
	ar0 = ar2;
	ar1 = 128;
L45:
	gr6 = gr0 << 2;
	gr3 = gr0 << 4;
	gr6 = gr6 + gr3;
	gr3 = gr6 << 6;
	gr6 = gr3 - gr6;
	gr6 = gr6 + gr0;
	gr3 = gr6 << 5;
	gr6 = gr6 + gr3;
	gr6 = gr6 << 3;
	gr0 = gr6 + gr0;
	gr3 = gr0 << 2;
	gr0 = gr0 + gr3;
	ar3 = gr0;
	ar3 = ar3 + 1013904223;
	gr0 = ar3;
	gr3 = gr0 >> 16;
	gr3 = gr3 and gr2;
	ar3 = gr3;
	ar3 = ar3 + -8;
	gr3 = ar3;
	gr7 = gr3;
	call ConvI32toF;
	gr3 = gr7;
	[ar0] = gr3;
	gr6 = gr0 << 2;
	gr3 = gr0 << 4;
	gr6 = gr6 + gr3;
	gr3 = gr6 << 6;
	gr6 = gr3 - gr6;
	gr6 = gr6 + gr0;
	gr3 = gr6 << 5;
	gr6 = gr6 + gr3;
	gr6 = gr6 << 3;
	gr0 = gr6 + gr0;
	gr3 = gr0 << 2;
	gr0 = gr0 + gr3;
	ar3 = gr0;
	ar3 = ar3 + 1013904223;
	gr0 = ar3;
	gr3 = gr0 >> 16;
	gr3 = gr3 and gr2;
	ar3 = gr3;
	ar3 = ar3 + -8;
	gr3 = ar3;
	gr7 = gr3;
	call ConvI32toF;
	gr3 = gr7;
	[ar0 + 1] = gr3;
	ar1 = ar1 + -1;
	ar0 = ar0 + 2;
	gr3 = ar1;
	gr3;
	if <>0 goto L45;
	ar2 = ar2 + 256;
	gr3 = ar2;
	gr3 - gr1;
	if =0 goto L46;
	ar0 = ar2;
	ar1 = 128;
	goto L45;
L56:
	ar1 = _template_buf;
	gr0 = 03f800000h;
	ar0 = ar1;
	ar5 = 128;
L42:
	[ar0] = gr0;
	[ar0 + 1] = gr0;
	ar5 = ar5 + -1;
	ar0 = ar0 + 2;
	gr2 = ar5;
	gr2;
	if <>0 goto L42;
	ar1 = ar1 + 256;
	gr3 = ar1;
	gr3 - gr1;
	if =0 goto L46;
	ar0 = ar1;
	ar5 = 128;
	goto L42;
L55:
	ar1 = _template_buf + 1;
	ar0 = ar1 + 32768;
	gr2 = ar0;
	gr1 = 03f800000h;
	gr0 = 128;
	ar5 = ar1;
	ar0 = 0;
L38:
	[ar5] = gr1;
	ar0 = ar0 + 1;
	ar5 = ar5 + 2;
	gr3 = ar0;
	gr3 - gr0;
	if <>0 goto L38;
	ar1 = ar1 + 256;
	gr3 = ar1;
	gr3 - gr2;
	if =0 goto L46;
	ar5 = ar1;
	ar0 = 0;
	goto L38;
L54:
	gr0 = 03f800000h;
	[_template_buf + 64] = gr0;
	gr7 = 0;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar6,gr6;
	return;
	.size	_fill_case.constprop.2, .-_fill_case.constprop.2
	.align 1;
LC0:
	.string	"Usage: solution.abs check <case_id>\n       solution.abs bench <case_id> <iters>\n"
	.align 1;
LC1:
	.string	"Unknown case %d\n"
	.align 1;
LC2:
	.string	"Mismatch at (%d,%d): got %d %d expected %d %d\n"
	.align 1;
LC3:
	.string	"OK case-%d\n"
	.align 1;
LC4:
	.string	"Usage: solution.abs bench <case_id> <iters>\n"
	.align 1;
LC5:
	.string	"Unknown mode: %s\n"
	.align 1;
.globl __main
	.type	__main, @function
__main:
	push ar6,gr6;
	push ar5,gr5;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	ar7 = ar7 + 40;
	ar6 = ar7 + -40;
	ar0 = ar6 + -16;
	gr0 = [ar0 + -1];
	ar1 = [ar0 + -2];
	gr1 = 2;
	gr0 - gr1;
	if v<= goto L171;
	ar0 = ar7 + -1;
	gr1 = [ar1 + 2];
	[ar0] = gr1;
	call _atoi;
	gr4 = gr7;
	ar2 = [ar1 + 1];
	gr1 = [ar2];
	gr2 = 99;
	gr1 - gr2;
	if =0 goto L172;
	gr2 = 98;
	gr1 - gr2;
	if =0 goto L173;
	ar1 = [__impure_ptr];
	ar1 = [ar1 + 3];
	[ar0] = ar1;
	gr0 = LC5;
	[ar7 + -2] = gr0;
	[ar7 + -3] = ar2;
	call _fprintf;
L165:
	gr2 = 1;
L62:
	gr7 = gr2;
	ar7 = ar7 + -40;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L173:
	gr1 = 3;
	gr0 - gr1;
	if =0 goto L174;
	ar1 = [ar1 + 3];
	[ar0] = ar1;
	call _atoi;
	gr3 = gr7;
	[ar0] = gr4;
	call _fill_case.constprop.2;
	gr7;
	if <>0 goto L165;
	gr3;
	if v<= goto L175;
	ar1 = 0;
	[ar6] = ar1;
	ar2 = ar1;
	ar5 = ar7 + -2;
	[ar6 + 2] = ar5;
	ar5 = ar7 + -3;
	[ar6 + 6] = ar5;
	gr1 = _matrix_buf;
	gr0 = _template_buf;
	gr2 = gr1;
L123:
	gr5 = 0;
	gr7 = 32768;
L129:
	gr4 = gr1 + gr5;
	ar3 = gr4;
	gr6 = gr0 + gr5;
	ar5 = gr6;
	gr6 = [ar5];
	[ar3] = gr6;
	ar5 = [ar5 + 1];
	[ar3 + 1] = ar5;
	ar3 = gr5;
	ar3 = ar3 + 2;
	gr5 = ar3;
	gr5 - gr7;
	if <>0 goto L129;
	[ar0] = gr2;
	gr5 = _twiddle_cos;
	ar5 = [ar6 + 2];
	[ar5] = gr5;
	gr5 = _twiddle_sin;
	ar4 = [ar6 + 6];
	[ar4] = gr5;
	call _contest_fft_2d;
	gr5 = 16383;
	gr6 = ar1;
	gr5 = gr6 and gr5;
	gr5 = gr5 + gr5;
	gr5 = gr2 + gr5;
	ar3 = gr5;
	gr5 = [ar3];
	gr4 = 000000000h;
	.align;
	push gr5;
	push gr4;
	call FCmp;
	ar7 -= 2;
	if v>= goto L151;
	gr4 = 03f000000h;
	gr7 = 80000000h;
	gr7 = gr7 xor gr4;
	delayed call FAdd;
	push gr5;
	push gr7;
	gr5 = gr7;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
	gr6 = [ar6];
	gr5 = gr6 + gr5;
	[ar6] = gr5;
	ar2 = ar2 + 1;
	ar1 = ar1 + 97;
	gr5 = ar2;
	gr5 - gr3;
	if v< goto L123;
L128:
	gr1 = [ar6];
	gr2 = 1;
	gr0 = 123456789;
	gr1 - gr0;
	if =0 goto L62;
	gr2 = 0;
	gr7 = gr2;
	ar7 = ar7 + -40;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L171:
	gr0 = LC0;
	[ar7 + -1] = gr0;
	gr2 = 1;
	[ar7 + -2] = gr2;
	gr0 = 80;
	[ar7 + -3] = gr0;
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar7 + -4] = ar0;
	call _fwrite;
L176:
	gr7 = gr2;
	ar7 = ar7 + -40;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L174:
	gr0 = LC4;
	[ar0] = gr0;
	gr2 = 1;
	[ar7 + -2] = gr2;
	gr0 = 44;
	[ar7 + -3] = gr0;
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar7 + -4] = ar0;
	call _fwrite;
	goto L176;
L151:
	gr4 = 03f000000h;
	delayed call FAdd;
	push gr5;
	push gr4;
	gr5 = gr7;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
	gr6 = [ar6];
	gr5 = gr6 + gr5;
	[ar6] = gr5;
	ar2 = ar2 + 1;
	ar1 = ar1 + 97;
	gr5 = ar2;
	gr5 - gr3;
	if v< goto L123;
	goto L128;
L172:
	[ar0] = gr7;
	call _fill_case.constprop.2;
	gr2 = gr7;
	gr7;
	if <>0 goto L177;
	gr3 = gr7;
	gr1 = _matrix_buf;
	gr0 = _template_buf;
	gr6 = 32768;
L67:
	ar3 = _matrix_buf;
	gr5 = gr1 + gr3;
	ar1 = gr5;
	gr5 = gr0 + gr3;
	ar2 = gr5;
	gr5 = [ar2];
	[ar1] = gr5;
	ar2 = [ar2 + 1];
	[ar1 + 1] = ar2;
	ar1 = gr3;
	ar1 = ar1 + 2;
	gr3 = ar1;
	gr3 - gr6;
	if <>0 goto L67;
	[ar0] = ar3;
	ar5 = ar7 + -2;
	[ar6 + 2] = ar5;
	gr3 = _twiddle_cos;
	[ar5] = gr3;
	ar4 = ar7 + -3;
	[ar6 + 6] = ar4;
	gr3 = _twiddle_sin;
	[ar4] = gr3;
	call _contest_fft_2d;
	gr3 = 5;
	gr4 - gr3;
	if =0 goto L68;
	ar4 = ar3 + 1;
	[ar6 + 10] = ar4;
	ar3 = 0;
	[ar6 + 14] = gr2;
	gr2 = ar3;
L69:
	ar2 = ar6 + 10;
	ar1 = [ar2];
	gr3 = 0;
	gr0 = gr2 and gr4;
	ar3 = gr0;
	[ar6 + 22] = ar0;
	ar0 = gr2;
L118:
	gr1 = [ar1 + -1];
	[ar6 + 26] = gr1;
	gr1 = [ar1];
	[ar6 + 18] = gr1;
	gr4;
	if =0 goto L131;
	gr1 = 1;
	gr4 - gr1;
	if =0 goto L178;
	gr6 = 2;
	gr4 - gr6;
	if =0 goto L179;
	gr0 = 3;
	gr4 - gr0;
	if =0 goto L180;
	gr7 = 4;
	gr4 - gr7;
	if =0 goto L181;
	gr0 = 6;
	gr4 - gr0;
	if =0 goto L182;
	gr0 = 7;
	gr4 - gr0;
	if =0 goto L183;
L141:
	gr0 = 0;
L163:
	gr1 = gr0;
L96:
	gr2 = [ar6 + 26];
	gr5 = 000000000h;
	.align;
	push gr2;
	push gr5;
	call FCmp;
	ar7 -= 2;
	if v>= goto L149;
	gr5 = [ar6 + 26];
	gr2 = 03f000000h;
	gr7 = 80000000h;
	gr7 = gr7 xor gr2;
	delayed call FAdd;
	push gr5;
	push gr7;
	gr5 = gr7;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
	[ar6 + 26] = gr5;
L110:
	gr2 = [ar6 + 18];
	gr5 = 000000000h;
	.align;
	push gr2;
	push gr5;
	call FCmp;
	ar7 -= 2;
	if v>= goto L150;
	gr5 = [ar6 + 18];
	gr2 = 03f000000h;
	gr7 = 80000000h;
	gr7 = gr7 xor gr2;
	delayed call FAdd;
	push gr5;
	push gr7;
	gr5 = gr7;
	[ar6 + 18] = gr5;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
L113:
	gr2 = [ar6 + 26];
	gr2 - gr1;
	if <>0 goto L162;
	gr5 - gr0;
	if <>0 goto L184;
	ar2 = gr3;
	ar2 = ar2 + 1;
	gr3 = ar2;
	ar1 = ar1 + 2;
	gr0 = 128;
	gr3 - gr0;
	if <>0 goto L118;
	gr2 = ar0;
	ar3 = ar6 + 22;
	ar0 = [ar3];
	ar4 = gr2;
	ar4 = ar4 + 1;
	gr2 = ar4;
	ar5 = [ar6 + 10];
	ar2 = ar5 + 256;
	[ar6 + 10] = ar2;
	gr2 - gr3;
	if <>0 goto L69;
	ar1 = ar6 + 14;
	gr2 = [ar1];
	gr0 = LC3;
	[ar0] = gr0;
	ar3 = [ar6 + 2];
	[ar3] = gr4;
	call _printf;
	gr7 = gr2;
	ar7 = ar7 + -40;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L183:
	gr2 = ar0;
	gr2;
	if <>0 goto L141;
	gr3;
	if <>0 goto L142;
	gr0 = 16384;
	gr1 = gr0;
	goto L96;
L175:
	gr3 = 1;
	ar1 = 0;
	[ar6] = ar1;
	ar2 = ar1;
	ar5 = ar7 + -2;
	[ar6 + 2] = ar5;
	ar5 = ar7 + -3;
	[ar6 + 6] = ar5;
	gr1 = _matrix_buf;
	gr0 = _template_buf;
	gr2 = gr1;
	goto L123;
L68:
	gr4 = 0;
	gr3 = _reference_buf;
	gr6 = 32768;
L71:
	gr5 = gr3 + gr4;
	ar1 = gr5;
	gr5 = gr0 + gr4;
	ar2 = gr5;
	gr5 = [ar2];
	[ar1] = gr5;
	ar2 = [ar2 + 1];
	[ar1 + 1] = ar2;
	ar1 = gr4;
	ar1 = ar1 + 2;
	gr4 = ar1;
	gr4 - gr6;
	if <>0 goto L71;
	gr0 = _reference_buf + 32768;
	ar1 = _reference_buf;
L74:
	[ar0] = ar1;
	call _fft1d_reference;
	ar1 = ar1 + 256;
	gr4 = ar1;
	gr4 - gr0;
	if <>0 goto L74;
	gr4 = _reference_buf + 256;
	ar2 = _reference_buf;
	gr5 = _vec_buf0 + 256;
	gr0 = _vec_buf0;
L73:
	ar1 = _vec_buf0;
	ar3 = ar2;
L79:
	gr6 = [ar3];
	[ar1] = gr6;
	gr6 = [ar3 + 1];
	[ar1 + 1] = gr6;
	ar3 = ar3 + 256;
	ar1 = ar1 + 2;
	gr6 = ar1;
	gr6 - gr5;
	if <>0 goto L79;
	[ar0] = gr0;
	call _fft1d_reference;
	ar3 = ar2;
	ar1 = _vec_buf0;
L77:
	gr6 = [ar1];
	[ar3] = gr6;
	gr6 = [ar1 + 1];
	[ar3 + 1] = gr6;
	ar1 = ar1 + 2;
	ar3 = ar3 + 256;
	gr6 = ar1;
	gr6 - gr5;
	if <>0 goto L77;
	ar2 = ar2 + 2;
	gr6 = ar2;
	gr6 - gr4;
	if <>0 goto L73;
	ar3 = 0;
	gr4 = ar3;
	gr0 = 03f000000h;
	ar1 = ar0;
	[ar6 + 14] = gr2;
	[ar6 + 18] = gr1;
L95:
	gr2 = [ar6 + 18];
	gr5 = ar3;
	gr1 = gr2 + gr5;
	ar0 = gr1;
	gr5 = [ar0];
	gr2 = 000000000h;
	.align;
	push gr5;
	push gr2;
	call FCmp;
	ar7 -= 2;
	if v>= goto L145;
	gr7 = 80000000h;
	gr7 = gr7 xor gr0;
	delayed call FAdd;
	push gr5;
	push gr7;
	gr5 = gr7;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
	ar2 = gr5;
L82:
	ar0 = gr1;
	gr1 = [ar0 + 1];
	gr2 = 000000000h;
	.align;
	push gr1;
	push gr2;
	call FCmp;
	ar7 -= 2;
	if v>= goto L146;
	gr7 = 80000000h;
	gr7 = gr7 xor gr0;
	delayed call FAdd;
	push gr1;
	push gr7;
	gr1 = gr7;
	gr7 = gr1;
	call ConvFtoI32;
	gr5 = gr7;
	[ar6 + 10] = gr5;
L85:
	gr6 = ar3;
	gr1 = gr3 + gr6;
	ar0 = gr1;
	gr5 = [ar0];
	gr2 = 000000000h;
	.align;
	push gr5;
	push gr2;
	call FCmp;
	ar7 -= 2;
	if v>= goto L147;
	gr7 = 80000000h;
	gr7 = gr7 xor gr0;
	delayed call FAdd;
	push gr5;
	push gr7;
	gr5 = gr7;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
L88:
	ar0 = gr1;
	gr1 = [ar0 + 1];
	gr2 = 000000000h;
	.align;
	push gr1;
	push gr2;
	call FCmp;
	ar7 -= 2;
	if v>= goto L148;
	gr7 = 80000000h;
	gr7 = gr7 xor gr0;
	delayed call FAdd;
	push gr1;
	push gr7;
	gr1 = gr7;
	gr7 = gr1;
	call ConvFtoI32;
	gr1 = gr7;
L91:
	gr6 = ar2;
	gr6 - gr5;
	if <>0 goto L161;
	gr2 = [ar6 + 10];
	gr2 - gr1;
	if <>0 goto L185;
	ar2 = gr4;
	ar2 = ar2 + 1;
	gr4 = ar2;
	ar3 = ar3 + 2;
	gr1 = 16384;
	gr4 - gr1;
	if <>0 goto L95;
	ar0 = ar6 + 14;
	gr2 = [ar0];
	gr0 = LC3;
	[ar1] = gr0;
	gr0 = 5;
	ar2 = [ar6 + 2];
	[ar2] = gr0;
	call _printf;
	gr7 = gr2;
	ar7 = ar7 + -40;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L150:
	gr5 = [ar6 + 18];
	gr2 = 03f000000h;
	delayed call FAdd;
	push gr5;
	push gr2;
	gr5 = gr7;
	[ar6 + 18] = gr5;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
	goto L113;
L149:
	gr5 = [ar6 + 26];
	gr2 = 03f000000h;
	delayed call FAdd;
	push gr5;
	push gr2;
	gr5 = gr7;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
	[ar6 + 26] = gr5;
	goto L110;
L131:
	gr0 = gr4;
	gr1 = 1;
	goto L96;
L148:
	delayed call FAdd;
	push gr1;
	push gr0;
	gr1 = gr7;
	gr7 = gr1;
	call ConvFtoI32;
	gr1 = gr7;
	goto L91;
L147:
	delayed call FAdd;
	push gr5;
	push gr0;
	gr5 = gr7;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
	goto L88;
L146:
	delayed call FAdd;
	push gr1;
	push gr0;
	gr1 = gr7;
	gr7 = gr1;
	call ConvFtoI32;
	gr5 = gr7;
	[ar6 + 10] = gr5;
	goto L85;
L145:
	delayed call FAdd;
	push gr5;
	push gr0;
	gr5 = gr7;
	gr7 = gr5;
	call ConvFtoI32;
	gr5 = gr7;
	ar2 = gr5;
	goto L82;
L178:
	gr2 = ar0;
	gr2;
	if <>0 goto L141;
	gr3;
	if <>0 goto L142;
	gr0 = gr3;
	gr1 = 16384;
	goto L96;
L179:
	gr0 = 64;
	gr5 = ar0;
	gr5 - gr0;
	if <>0 goto L141;
	gr0 = 0;
	gr3 - gr5;
	if <>0 goto L163;
	gr1 = 16384;
	goto L96;
L142:
	gr0 = ar0;
	gr1 = ar0;
	goto L96;
L180:
	gr0 = ar3;
	gr0 - gr1;
	if =0 goto L137;
L169:
	gr0 - gr6;
	if =0 goto L105;
	gr0;
	if =0 goto L96;
	gr0 = gr1;
	gr1 = 0;
	goto L96;
L181:
	gr0 = gr3 and gr0;
	gr0 - gr1;
	if <>0 goto L169;
L137:
	gr0 = -1;
	gr1 = 0;
	goto L96;
L105:
	gr0 = 0;
	gr1 = -1;
	goto L96;
L182:
	gr1 = ar0;
	gr1;
	if <>0 goto L141;
	gr3;
	if <>0 goto L142;
	gr0 = 16384;
	gr1 = gr3;
	goto L96;
L177:
	gr0 = LC1;
	[ar0] = gr0;
	[ar7 + -2] = gr4;
	call _printf;
	gr2 = 1;
	gr7 = gr2;
	ar7 = ar7 + -40;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L184:
	ar3 = ar0;
	ar1 = ar6 + 22;
	ar0 = [ar1];
	ar2 = gr2;
	gr1 = gr2;
	gr2 = LC2;
	[ar0] = gr2;
	ar5 = [ar6 + 2];
	[ar5] = ar3;
	ar1 = [ar6 + 6];
	[ar1] = gr3;
	[ar7 + -4] = ar2;
	[ar7 + -5] = gr5;
	[ar7 + -6] = gr1;
	[ar7 + -7] = gr0;
	call _printf;
	gr2 = 1;
L186:
	gr7 = gr2;
	ar7 = ar7 + -40;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L162:
	ar3 = ar0;
	ar5 = ar6 + 22;
	ar0 = [ar5];
	ar1 = ar6 + 26;
	ar2 = [ar1];
	gr2 = LC2;
	[ar0] = gr2;
	ar5 = [ar6 + 2];
	[ar5] = ar3;
	ar1 = [ar6 + 6];
	[ar1] = gr3;
	[ar7 + -4] = ar2;
	[ar7 + -5] = gr5;
	[ar7 + -6] = gr1;
	[ar7 + -7] = gr0;
	call _printf;
	gr2 = 1;
	goto L186;
L185:
	ar0 = ar1;
	ar1 = gr1;
	gr5 = ar2;
	gr0 = LC2;
	[ar0] = gr0;
	gr0 = gr4 A>> 7;
	ar3 = [ar6 + 2];
	[ar3] = gr0;
	gr0 = 127;
	gr4 = gr4 and gr0;
	ar3 = [ar6 + 6];
	[ar3] = gr4;
	[ar7 + -4] = ar2;
	ar0 = [ar6 + 10];
	[ar7 + -5] = ar0;
	[ar7 + -6] = gr5;
	[ar7 + -7] = ar1;
	call _printf;
	gr2 = 1;
L187:
	gr7 = gr2;
	ar7 = ar7 + -40;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L161:
	ar0 = ar1;
	ar1 = gr1;
	gr0 = LC2;
	[ar0] = gr0;
	gr0 = gr4 A>> 7;
	ar3 = [ar6 + 2];
	[ar3] = gr0;
	gr0 = 127;
	gr4 = gr4 and gr0;
	ar3 = [ar6 + 6];
	[ar3] = gr4;
	[ar7 + -4] = ar2;
	ar0 = [ar6 + 10];
	[ar7 + -5] = ar0;
	[ar7 + -6] = gr5;
	[ar7 + -7] = ar1;
	call _printf;
	gr2 = 1;
	goto L187;
	.size	__main, .-__main
	.data
	.align 3;
	.type	_twiddle_sin, @object
	.size	_twiddle_sin, 256
_twiddle_sin:
	.long	0
	.long	1028193071
	.long	1036565814
	.long	1041645699
	.long	1044891074
	.long	1048104908
	.long	1049927729
	.long	1051491540
	.long	1053028117
	.long	1054533760
	.long	1056004842
	.long	1057201213
	.long	1057896922
	.long	1058570176
	.long	1059219353
	.long	1059842890
	.long	1060439283
	.long	1061007097
	.long	1061544963
	.long	1062051586
	.long	1062525745
	.long	1062966298
	.long	1063372184
	.long	1063742424
	.long	1064076126
	.long	1064372488
	.long	1064630795
	.long	1064850424
	.long	1065030846
	.long	1065171628
	.long	1065272429
	.long	1065333007
	.long	1065353216
	.long	1065333007
	.long	1065272429
	.long	1065171628
	.long	1065030846
	.long	1064850424
	.long	1064630795
	.long	1064372488
	.long	1064076126
	.long	1063742424
	.long	1063372184
	.long	1062966298
	.long	1062525745
	.long	1062051586
	.long	1061544963
	.long	1061007097
	.long	1060439283
	.long	1059842890
	.long	1059219353
	.long	1058570176
	.long	1057896922
	.long	1057201213
	.long	1056004842
	.long	1054533760
	.long	1053028117
	.long	1051491540
	.long	1049927729
	.long	1048104908
	.long	1044891074
	.long	1041645699
	.long	1036565814
	.long	1028193071
	.align 3;
	.type	_twiddle_cos, @object
	.size	_twiddle_cos, 256
_twiddle_cos:
	.long	1065353216
	.long	1065333007
	.long	1065272429
	.long	1065171628
	.long	1065030846
	.long	1064850424
	.long	1064630795
	.long	1064372488
	.long	1064076126
	.long	1063742424
	.long	1063372184
	.long	1062966298
	.long	1062525745
	.long	1062051586
	.long	1061544963
	.long	1061007097
	.long	1060439283
	.long	1059842890
	.long	1059219353
	.long	1058570176
	.long	1057896922
	.long	1057201213
	.long	1056004842
	.long	1054533760
	.long	1053028117
	.long	1051491540
	.long	1049927729
	.long	1048104908
	.long	1044891074
	.long	1041645699
	.long	1036565814
	.long	1028193071
	.long	0
	.long	3175676719
	.long	3184049462
	.long	3189129347
	.long	3192374722
	.long	3195588556
	.long	3197411377
	.long	3198975188
	.long	3200511765
	.long	3202017408
	.long	3203488490
	.long	3204684861
	.long	3205380570
	.long	3206053824
	.long	3206703001
	.long	3207326538
	.long	3207922931
	.long	3208490745
	.long	3209028611
	.long	3209535234
	.long	3210009393
	.long	3210449946
	.long	3210855832
	.long	3211226072
	.long	3211559774
	.long	3211856136
	.long	3212114443
	.long	3212334072
	.long	3212514494
	.long	3212655276
	.long	3212756077
	.long	3212816655
	.section	.bss.emi,"aw",@nobits
	.align 3;
	.type	_vec_buf1, @object
	.size	_vec_buf1, 1024
_vec_buf1:
	.zero	1024
	.align 3;
	.type	_vec_buf0, @object
	.size	_vec_buf0, 1024
_vec_buf0:
	.zero	1024
	.align 3;
	.type	_reference_buf, @object
	.size	_reference_buf, 131072
_reference_buf:
	.zero	131072
	.align 3;
	.type	_template_buf, @object
	.size	_template_buf, 131072
_template_buf:
	.zero	131072
	.align 3;
	.type	_matrix_buf, @object
	.size	_matrix_buf, 131072
_matrix_buf:
	.zero	131072
	.ident	"GCC: (NMC SDK GCC ) 4.8.3"
