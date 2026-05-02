	.file	"harness_main.c"
	.text
	.align 1;
LC0:
	.string	"Usage: nmc-qemu solution.abs <base> <exp> <mod> [iterations]\n"
	.align 1;
LC1:
	.string	"Invalid base: %s\n"
	.align 1;
LC2:
	.string	"Invalid exp: %s\n"
	.align 1;
LC3:
	.string	"Invalid mod: %s\n"
	.align 1;
LC4:
	.string	"%u\n"
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
	ar7 = ar7 + 32;
	ar6 = ar7 + -32;
	ar3 = ar6 + -16;
	gr4 = [ar3 + -1];
	gr6 = [ar3 + -2];
	[ar6 + 22] = gr6;
	gr5 = 3;
	gr4 - gr5;
	if v<= goto L15;
	ar1 = [ar6 + 22];
	ar0 = ar1 + 1;
	gr5 = ar0;
	ar0 = ar7 + -1;
	ar1 = gr5;
	gr6 = [ar1];
	[ar0] = gr6;
	ar1 = ar7 + -2;
	[ar1] = ar6;
	ar2 = ar7 + -3;
	gr6 = 0;
	[ar2] = gr6;
	call _strtoul;
	gr2 = gr7;
	ar5 = [ar6];
	gr6 = [ar5];
	gr6;
	if <>0 goto L16;
	ar5 = [ar6 + 22];
	ar4 = ar5 + 2;
	gr7 = [ar4];
	[ar0] = gr7;
	[ar1] = ar6;
	[ar2] = gr6;
	call _strtoul;
	gr0 = gr7;
	ar5 = [ar6];
	gr6 = [ar5];
	gr6;
	if <>0 goto L17;
	ar3 = [ar6 + 22];
	ar5 = ar3 + 3;
	[ar6 + 18] = ar5;
	gr5 = [ar5];
	[ar0] = gr5;
	[ar1] = ar6;
	[ar2] = gr6;
	call _strtoul;
	[ar6 + 10] = gr7;
	[ar6 + 11] = gr6;
	ar3 = ar6 + 10;
	gr5 = [ar3];
	ar5 = [ar6];
	gr6 = [ar5];
	gr6;
	if <>0 goto L18;
	gr6 = 4;
	gr4 - gr6;
	if <>0 goto L19;
	gr4 = 1;
L8:
	gr6 = 0;
	[ar6 + 8] = gr6;
	ar3 = gr6;
L10:
	[ar0] = gr2;
	[ar1] = gr0;
	[ar2] = gr5;
	call _contest_pow_mod;
	[ar6 + 8] = gr7;
	ar3 = ar3 + 1;
	gr6 = ar3;
	gr4 - gr6;
	if v> goto L10;
	gr0 = [ar6 + 8];
	gr1 = LC4;
	[ar0] = gr1;
	[ar1] = gr0;
	call _printf;
	gr7 = 0;
	ar7 = ar7 + -32;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L16:
	ar3 = [__impure_ptr];
	ar3 = [ar3 + 3];
	[ar0] = ar3;
	gr0 = LC1;
	[ar1] = gr0;
	ar3 = gr5;
	gr0 = [ar3];
	[ar2] = gr0;
	call _fprintf;
	gr7 = 1;
	ar7 = ar7 + -32;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L15:
	gr0 = LC0;
	[ar7 + -1] = gr0;
	gr0 = 1;
	[ar7 + -2] = gr0;
	gr1 = 61;
	[ar7 + -3] = gr1;
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar7 + -4] = ar0;
	call _fwrite;
	gr7 = gr0;
	ar7 = ar7 + -32;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L19:
	ar4 = [ar6 + 22];
	ar3 = [ar4 + 4];
	[ar0] = ar3;
	call _atoi;
	gr4 = gr7;
	gr7;
	if v> goto L8;
	gr4 = 1;
	goto L8;
L17:
	ar3 = [__impure_ptr];
	ar3 = [ar3 + 3];
	[ar0] = ar3;
	gr0 = LC2;
	[ar1] = gr0;
	gr0 = [ar4];
	[ar2] = gr0;
	call _fprintf;
	gr7 = 1;
	ar7 = ar7 + -32;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L18:
	ar3 = [__impure_ptr];
	ar3 = [ar3 + 3];
	[ar0] = ar3;
	gr0 = LC3;
	[ar1] = gr0;
	ar4 = [ar6 + 18];
	gr0 = [ar4];
	[ar2] = gr0;
	call _fprintf;
	gr7 = 1;
	ar7 = ar7 + -32;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size	__main, .-__main
	.ident	"GCC: (NMC SDK GCC ) 4.8.3"
