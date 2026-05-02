	.file	"harness_main.c"
	.text
	.align 1;
LC0:
	.string	"Usage: nmc-qemu solution.abs <n> [iterations]\n"
	.align 1;
LC1:
	.string	"Invalid n: %s\n"
	.align 1;
LC2:
	.string	"%d\n"
	.align 1;
.globl __main
	.type	__main, @function
__main:
	push ar6,gr6;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	ar7 = ar7 + 16;
	ar6 = ar7 + -16;
	ar0 = ar6 + -14;
	gr1 = [ar0 + -1];
	ar0 = [ar0 + -2];
	gr0 = 1;
	gr1 - gr0;
	if v<= goto L14;
	ar1 = ar0 + 1;
	gr4 = ar1;
	ar1 = ar7 + -1;
	ar2 = gr4;
	gr6 = [ar2];
	[ar1] = gr6;
	ar3 = ar7 + -2;
	[ar3] = ar6;
	ar2 = ar7 + -3;
	gr6 = 10;
	[ar2] = gr6;
	call _strtol;
	gr2 = gr7;
	ar5 = [ar6];
	gr3 = [ar5];
	gr3;
	if <>0 goto L4;
	gr7;
	if v< goto L4;
	gr3 = 2;
	gr1 - gr3;
	if <>0 goto L15;
	gr3 = gr0;
L7:
	ar2 = ar6 + 8;
	ar0 = 0;
	[ar2] = ar0;
L9:
	[ar1] = gr2;
	call _contest_last_nonzero_digit_factorial;
	[ar2] = gr7;
	ar0 = ar0 + 1;
	gr0 = ar0;
	gr3 - gr0;
	if v> goto L9;
	gr0 = [ar2];
	gr1 = LC2;
	[ar1] = gr1;
	[ar3] = gr0;
	call _printf;
	gr7 = 0;
	ar7 = ar7 + -16;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar6,gr6;
	return;
L4:
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar1] = ar0;
	gr0 = LC1;
	[ar3] = gr0;
	ar0 = gr4;
	gr0 = [ar0];
	[ar2] = gr0;
	call _fprintf;
	gr7 = 1;
	ar7 = ar7 + -16;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar6,gr6;
	return;
L15:
	ar0 = [ar0 + 2];
	[ar1] = ar0;
	call _atoi;
	gr3 = gr7;
	gr7;
	if v> goto L7;
	gr3 = gr0;
	goto L7;
L14:
	gr1 = LC0;
	[ar7 + -1] = gr1;
	[ar7 + -2] = gr0;
	gr1 = 46;
	[ar7 + -3] = gr1;
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar7 + -4] = ar0;
	call _fwrite;
	gr7 = gr0;
	ar7 = ar7 + -16;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar6,gr6;
	return;
	.size	__main, .-__main
	.ident	"GCC: (NMC SDK GCC ) 4.8.3"
