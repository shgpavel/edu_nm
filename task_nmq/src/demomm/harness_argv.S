	.file	"harness_argv.c"
	.text
	.align 1;
LC0:
	.string	"Usage: solution.abs gen m k n iters\n"
	.align 1;
LC1:
	.string	"Usage: solution.abs m k n a00 ... b00 ...\n       solution.abs gen m k n iters\n"
	.align 1;
LC2:
	.string	"Need %d args, got %d\n"
	.align 1;
LC3:
	.string	"%d"
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
	ar7 = ar7 + 28;
	ar6 = ar7 + -28;
	ar2 = ar6 + -16;
	gr2 = [ar2 + -1];
	ar2 = [ar2 + -2];
	gr0 = 1;
	gr2 - gr0;
	if v<= goto L2;
	ar1 = [ar2 + 1];
	gr1 = [ar1];
	gr3 = 103;
	gr1 - gr3;
	if =0 goto L42;
L3:
	gr0 = 3;
	gr2 - gr0;
	if v<= goto L2;
	ar0 = ar7 + -1;
	[ar0] = ar1;
	call _atoi;
	gr0 = ar6;
	[gr0] = gr7;
	gr0 = [ar2 + 2];
	[ar0] = gr0;
	call _atoi;
	[ar6 + 4] = gr7;
	gr0 = [ar2 + 3];
	[ar0] = gr0;
	call _atoi;
	gr0 = gr7;
	gr3 = [ar6];
	gr1 = [ar6 + 4];
	delayed call Mul32;
	push gr3;
	push gr1;
	gr4 = gr7;
	ar1 = gr4;
	ar4 = ar1 + 4;
	[ar6 + 8] = ar4;
	delayed call Mul32;
	push gr1;
	push gr0;
	gr5 = gr7;
	gr1 = ar4;
	gr3 = gr1 + gr5;
	gr2 - gr3;
	if v< goto L15;
	gr4;
	if v<= goto L19;
	ar1 = ar2 + 4;
	ar3 = _A;
	gr6 = ar1;
	gr1 = gr6 + gr4;
L20:
	gr6 = [ar1];
	[ar0] = gr6;
	call _atoi;
	[ar3] = gr7;
	ar1 = ar1 + 1;
	ar3 = ar3 + 1;
	gr7 = ar1;
	gr7 - gr1;
	if <>0 goto L20;
L19:
	gr5;
	if v<= goto L18;
	gr1 = ar2;
	gr6 = [ar6 + 8];
	gr1 = gr1 + gr6;
	ar3 = gr1;
	ar1 = _B;
	ar4 = ar2 + 4;
	gr1 = ar4;
	gr4 = gr1 + gr4;
	gr5 = gr4 + gr5;
L23:
	gr1 = [ar3];
	[ar0] = gr1;
	call _atoi;
	[ar1] = gr7;
	ar3 = ar3 + 1;
	ar1 = ar1 + 1;
	gr1 = ar3;
	gr1 - gr5;
	if <>0 goto L23;
L18:
	gr2 - gr3;
	if v> goto L21;
	gr2 = 1;
	[ar6 + 16] = gr2;
	gr1 = gr2;
	gr5 = 0;
	ar3 = ar7 + -2;
	ar4 = ar7 + -3;
	[ar6 + 12] = ar4;
	ar4 = ar7 + -4;
	[ar6 + 8] = ar4;
	gr4 = _A;
	gr3 = _B;
	gr2 = _C;
	ar2 = ar7 + -5;
	ar1 = ar7 + -6;
	gr6 = ar0;
	ar0 = gr5;
	gr5 = gr6;
L24:
	ar4 = gr5;
	[ar4] = gr4;
	[ar3] = gr3;
	ar4 = [ar6 + 12];
	[ar4] = gr2;
	gr7 = [ar6];
	ar5 = [ar6 + 8];
	[ar5] = gr7;
	ar4 = [ar6 + 4];
	[ar2] = ar4;
	[ar1] = gr0;
	call _contest_matmul_i32;
	ar0 = ar0 + 1;
	gr6 = ar0;
	gr1 - gr6;
	if v> goto L24;
	ar0 = gr5;
	gr1 = [ar6 + 16];
	gr1;
	if =0 goto L13;
	gr2 = [ar6];
	gr2;
	if v<= goto L13;
	gr2 = _C;
	ar1 = 0;
	gr3 = LC3;
	gr4 = 32;
	gr5 = 10;
	ar2 = [ar6];
L26:
	gr0;
	if v<= goto L29;
	gr1 = 0;
L30:
	[ar0] = gr3;
	gr6 = gr2 + gr1;
	ar5 = gr6;
	gr6 = [ar5];
	[ar3] = gr6;
	call _printf;
	ar4 = gr1;
	ar4 = ar4 + 1;
	gr1 = ar4;
	gr1 - gr0;
	if =0 goto L29;
L28:
	gr1;
	if =0 goto L30;
	[ar0] = gr4;
	call _putchar;
	[ar0] = gr3;
	gr6 = gr2 + gr1;
	ar5 = gr6;
	gr6 = [ar5];
	[ar3] = gr6;
	call _printf;
	ar4 = gr1;
	ar4 = ar4 + 1;
	gr1 = ar4;
	gr1 - gr0;
	if <>0 goto L28;
L29:
	[ar0] = gr5;
	call _putchar;
	ar1 = ar1 + 1;
	gr2 = gr2 + gr0;
	gr1 = ar1;
	gr6 = ar2;
	gr1 - gr6;
	if <>0 goto L26;
L13:
	gr7 = 0;
	ar7 = ar7 + -28;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L42:
	gr3 = [ar1 + 1];
	gr1 = 101;
	gr3 - gr1;
	if <>0 goto L3;
	gr3 = [ar1 + 2];
	gr1 = 110;
	gr3 - gr1;
	if <>0 goto L3;
	gr1 = 5;
	gr2 - gr1;
	if v<= goto L43;
	ar0 = ar7 + -1;
	gr0 = [ar2 + 2];
	[ar0] = gr0;
	call _atoi;
	gr3 = gr7;
	gr0 = [ar2 + 3];
	[ar0] = gr0;
	call _atoi;
	gr2 = gr7;
	gr0 = [ar2 + 4];
	[ar0] = gr0;
	call _atoi;
	gr4 = gr7;
	ar2 = [ar2 + 5];
	[ar0] = ar2;
	call _atoi;
	gr5 = gr7;
	delayed call Mul32;
	push gr3;
	push gr2;
	gr0 = gr7;
	gr0;
	if v<= goto L31;
	ar1 = _A;
	gr1 = ar1;
	gr6 = gr0 + gr1;
	gr0 = 42;
L7:
	gr7 = gr0 << 2;
	gr1 = gr0 << 4;
	gr1 = gr7 + gr1;
	gr7 = gr1 << 6;
	gr1 = gr7 - gr1;
	gr1 = gr1 + gr0;
	gr7 = gr1 << 5;
	gr1 = gr1 + gr7;
	gr1 = gr1 << 3;
	gr0 = gr1 + gr0;
	gr1 = gr0 << 2;
	gr0 = gr0 + gr1;
	ar2 = gr0;
	ar2 = ar2 + 1013904223;
	gr0 = ar2;
	gr1 = gr0 A>> 16;
	[ar1] = gr1;
	ar1 = ar1 + 1;
	gr1 = ar1;
	gr1 - gr6;
	if =0 goto L6;
	gr7 = gr0 << 2;
	gr1 = gr0 << 4;
	gr1 = gr7 + gr1;
	gr7 = gr1 << 6;
	gr1 = gr7 - gr1;
	gr1 = gr1 + gr0;
	gr7 = gr1 << 5;
	gr1 = gr1 + gr7;
	gr1 = gr1 << 3;
	gr0 = gr1 + gr0;
	gr1 = gr0 << 2;
	gr0 = gr0 + gr1;
	ar2 = gr0;
	ar2 = ar2 + 1013904223;
	gr0 = ar2;
	gr1 = gr0 A>> 16;
	[ar1] = gr1;
	ar1 = ar1 + 1;
	gr1 = ar1;
	gr1 - gr6;
	if <>0 goto L7;
	goto L6;
L2:
	gr0 = LC1;
	[ar7 + -1] = gr0;
	gr0 = 1;
	[ar7 + -2] = gr0;
	gr1 = 78;
	[ar7 + -3] = gr1;
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar7 + -4] = ar0;
	call _fwrite;
	gr7 = gr0;
L44:
	ar7 = ar7 + -28;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L15:
	ar1 = [__impure_ptr];
	ar1 = [ar1 + 3];
	[ar0] = ar1;
	gr0 = LC2;
	[ar7 + -2] = gr0;
	[ar7 + -3] = gr3;
	[ar7 + -4] = gr2;
	call _fprintf;
	gr7 = 1;
	ar7 = ar7 + -28;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L21:
	gr6 = ar2;
	gr6 = gr6 + gr3;
	ar2 = gr6;
	gr1 = [ar2];
	[ar0] = gr1;
	call _atoi;
	gr1 = gr7;
	gr7;
	if v<= goto L13;
	gr3 = 0;
	[ar6 + 16] = gr3;
	gr5 = 0;
	ar3 = ar7 + -2;
	ar4 = ar7 + -3;
	[ar6 + 12] = ar4;
	ar4 = ar7 + -4;
	[ar6 + 8] = ar4;
	gr4 = _A;
	gr3 = _B;
	gr2 = _C;
	ar2 = ar7 + -5;
	ar1 = ar7 + -6;
	gr6 = ar0;
	ar0 = gr5;
	gr5 = gr6;
	goto L24;
L43:
	gr1 = LC0;
	[ar7 + -1] = gr1;
	[ar7 + -2] = gr0;
	gr1 = 36;
	[ar7 + -3] = gr1;
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar7 + -4] = ar0;
	call _fwrite;
	gr7 = gr0;
	goto L44;
L31:
	gr0 = 42;
L6:
	delayed call Mul32;
	push gr2;
	push gr4;
	gr1 = gr7;
	gr1;
	if v<= goto L10;
	ar1 = _B;
	gr7 = ar1;
	gr6 = gr1 + gr7;
L11:
	gr7 = gr0 << 2;
	gr1 = gr0 << 4;
	gr1 = gr7 + gr1;
	gr7 = gr1 << 6;
	gr1 = gr7 - gr1;
	gr1 = gr1 + gr0;
	gr7 = gr1 << 5;
	gr1 = gr1 + gr7;
	gr1 = gr1 << 3;
	gr0 = gr1 + gr0;
	gr1 = gr0 << 2;
	gr0 = gr0 + gr1;
	ar2 = gr0;
	ar2 = ar2 + 1013904223;
	gr0 = ar2;
	gr1 = gr0 A>> 16;
	[ar1] = gr1;
	ar1 = ar1 + 1;
	gr1 = ar1;
	gr1 - gr6;
	if <>0 goto L11;
L10:
	ar1 = 0;
	ar2 = ar7 + -2;
	gr1 = ar2;
	ar3 = ar7 + -3;
	gr0 = ar3;
	ar3 = ar7 + -4;
	ar2 = ar7 + -5;
L9:
	gr6 = ar1;
	gr6 - gr5;
	if v>= goto L13;
	gr6 = _A;
	[ar0] = gr6;
	gr6 = _B;
	ar4 = gr1;
	[ar4] = gr6;
	gr6 = _C;
	ar5 = gr0;
	[ar5] = gr6;
	[ar3] = gr3;
	[ar2] = gr2;
	[ar7 + -6] = gr4;
	call _contest_matmul_i32;
	ar1 = ar1 + 1;
	goto L9;
	.size	__main, .-__main
	.section	.bss.emi,"aw",@nobits
	.align 3;
	.type	_C, @object
	.size	_C, 4194304
_C:
	.zero	4194304
	.align 3;
	.type	_B, @object
	.size	_B, 4194304
_B:
	.zero	4194304
	.align 3;
	.type	_A, @object
	.size	_A, 4194304
_A:
	.zero	4194304
	.ident	"GCC: (NMC SDK GCC ) 4.8.3"
