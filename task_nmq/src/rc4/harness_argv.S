	.file	"harness_argv.c"
	.text
	.align 1;
LC0:
	.string	"Usage: solution.abs key_len k0..kN data_len d0..dM [iters]\n       solution.abs gen key_len data_len iters\n"
	.align 1;
LC1:
	.string	"Usage: solution.abs gen key_len data_len iters\n"
	.align 1;
LC2:
	.string	"%d "
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
	ar7 = ar7 + 276;
	ar6 = ar7 + -276;
	ar0 = ar6 + -16;
	gr1 = [ar0 + -1];
	[ar6 + 266] = gr1;
	gr4 = [ar0 + -2];
	gr0 = 1;
	gr1 - gr0;
	if v<= goto L34;
	ar0 = gr4;
	ar1 = [ar0 + 1];
	gr1 = [ar1];
	gr2 = 103;
	gr1 - gr2;
	if =0 goto L35;
L4:
	ar0 = ar7 + -1;
	[ar0] = ar1;
	call _atoi;
	[ar6 + 258] = gr7;
	gr7;
	if v<= goto L29;
	ar1 = gr7;
	ar1 = ar1 + 2;
	gr5 = ar1;
	gr0 = 2;
	gr1 = _key_buf;
L17:
	gr2 = gr4 + gr0;
	ar1 = gr2;
	gr2 = [ar1];
	[ar0] = gr2;
	call _atoi;
	gr2 = gr1 + gr0;
	ar1 = gr2;
	[ar1 + -2] = gr7;
	ar2 = gr0;
	ar2 = ar2 + 1;
	gr0 = ar2;
	gr0 - gr5;
	if <>0 goto L17;
L15:
	ar3 = gr0;
	ar5 = ar3 + 1;
	[ar6 + 262] = ar5;
	gr5 = gr4 + gr0;
	ar1 = gr5;
	gr0 = [ar1];
	[ar0] = gr0;
	call _atoi;
	gr2 = gr7;
	gr7;
	if v<= goto L18;
	ar2 = gr7;
	ar2 = ar2 + 1;
	gr3 = ar2;
	gr1 = 1;
	gr0 = _input_buf;
L20:
	gr6 = gr5 + gr1;
	ar2 = gr6;
	gr6 = [ar2];
	[ar0] = gr6;
	call _atoi;
	gr6 = gr0 + gr1;
	ar1 = gr6;
	[ar1 + -1] = gr7;
	ar2 = gr1;
	ar2 = ar2 + 1;
	gr1 = ar2;
	gr1 - gr3;
	if <>0 goto L20;
	gr1 = [ar6 + 262];
	gr1 = gr2 + gr1;
	[ar6 + 262] = gr1;
L18:
	gr0 = [ar6 + 266];
	gr1 = [ar6 + 262];
	gr0 - gr1;
	if v> goto L21;
	gr4 = 1;
	ar3 = 0;
	ar1 = ar7 + -2;
	ar2 = ar7 + -3;
	ar4 = ar7 + -4;
	gr3 = ar4;
	gr1 = _key_buf;
	gr0 = _input_buf;
	gr5 = _output_buf;
L22:
	[ar0] = ar6;
	[ar1] = gr1;
	ar4 = [ar6 + 258];
	[ar2] = ar4;
	call _contest_rc4_init;
	[ar0] = ar6;
	[ar1] = gr0;
	[ar2] = gr5;
	ar5 = gr3;
	[ar5] = gr2;
	call _contest_rc4_crypt;
	ar3 = ar3 + 1;
	gr6 = ar3;
	gr4 - gr6;
	if v> goto L22;
	gr0 = 1;
	gr4 - gr0;
	if v<= goto L24;
	gr0 = [ar6 + 266];
	gr1 = [ar6 + 262];
	gr0 - gr1;
	if v> goto L30;
L24:
	gr2;
	if v<= goto L27;
	ar2 = _output_buf;
	gr6 = ar2;
	gr2 = gr2 + gr6;
	ar1 = ar7 + -2;
	gr1 = LC2;
L28:
	[ar0] = gr1;
	gr0 = [ar2];
	[ar1] = gr0;
	call _printf;
	ar2 = ar2 + 1;
	gr0 = ar2;
	gr0 - gr2;
	if <>0 goto L28;
L27:
	gr0 = 10;
	[ar0] = gr0;
	call _putchar;
	gr7 = 0;
	ar7 = ar7 + -276;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L35:
	gr2 = [ar1 + 1];
	gr1 = 101;
	gr2 - gr1;
	if <>0 goto L4;
	gr2 = [ar1 + 2];
	gr1 = 110;
	gr2 - gr1;
	if <>0 goto L4;
	gr1 = 4;
	gr2 = [ar6 + 266];
	gr2 - gr1;
	if v> goto L5;
	gr1 = LC1;
	[ar7 + -1] = gr1;
	[ar7 + -2] = gr0;
	gr1 = 47;
	[ar7 + -3] = gr1;
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar7 + -4] = ar0;
	call _fwrite;
	gr7 = gr0;
L36:
	ar7 = ar7 + -276;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L30:
	gr7 = 0;
	ar7 = ar7 + -276;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
L21:
	gr0 = [ar6 + 262];
	gr4 = gr4 + gr0;
	ar1 = gr4;
	gr0 = [ar1];
	[ar0] = gr0;
	call _atoi;
	gr4 = gr7;
	gr7;
	if v<= goto L24;
	ar3 = 0;
	ar1 = ar7 + -2;
	ar2 = ar7 + -3;
	ar4 = ar7 + -4;
	gr3 = ar4;
	gr1 = _key_buf;
	gr0 = _input_buf;
	gr5 = _output_buf;
	goto L22;
L34:
	gr1 = LC0;
	[ar7 + -1] = gr1;
	[ar7 + -2] = gr0;
	gr1 = 106;
	[ar7 + -3] = gr1;
	ar0 = [__impure_ptr];
	ar0 = [ar0 + 3];
	[ar7 + -4] = ar0;
	call _fwrite;
	gr7 = gr0;
	goto L36;
L5:
	ar0 = ar7 + -1;
	ar1 = gr4;
	gr0 = [ar1 + 2];
	[ar0] = gr0;
	call _atoi;
	gr5 = gr7;
	gr0 = [ar1 + 3];
	[ar0] = gr0;
	call _atoi;
	ar1 = gr7;
	ar2 = gr4;
	gr4 = [ar2 + 4];
	[ar0] = gr4;
	call _atoi;
	gr4 = gr7;
	gr5;
	if v<= goto L10;
	gr0 = 0;
	gr1 = _key_buf;
L11:
	gr2 = gr1 + gr0;
	ar2 = gr2;
	[ar2] = gr0;
	ar3 = gr0;
	ar3 = ar3 + 1;
	gr0 = ar3;
	gr0 - gr5;
	if <>0 goto L11;
L10:
	gr0 = ar1;
	gr0;
	if v<= goto L9;
	gr1 = 0;
	gr0 = _input_buf;
	gr6 = 255;
L13:
	gr2 = gr0 + gr1;
	gr3 = gr1 and gr6;
	ar2 = gr2;
	[ar2] = gr3;
	ar3 = gr1;
	ar3 = ar3 + 1;
	gr1 = ar3;
	gr2 = ar1;
	gr1 - gr2;
	if <>0 goto L13;
L9:
	ar2 = 0;
	ar4 = ar7 + -2;
	gr0 = ar4;
	gr3 = _key_buf;
	ar3 = ar7 + -3;
	ar5 = ar7 + -4;
	gr2 = ar5;
L12:
	gr1 = ar2;
	gr1 - gr4;
	if v>= goto L30;
	[ar0] = ar6;
	ar4 = gr0;
	[ar4] = gr3;
	[ar3] = gr5;
	call _contest_rc4_init;
	[ar0] = ar6;
	gr1 = _input_buf;
	[ar4] = gr1;
	gr1 = _output_buf;
	[ar3] = gr1;
	ar5 = gr2;
	[ar5] = ar1;
	call _contest_rc4_crypt;
	ar2 = ar2 + 1;
	goto L12;
L29:
	gr0 = 2;
	goto L15;
	.size	__main, .-__main
	.local	_output_buf
	.comm	_output_buf,16384,8
	.local	_input_buf
	.comm	_input_buf,16384,8
	.local	_key_buf
	.comm	_key_buf,1024,8
	.ident	"GCC: (NMC SDK GCC ) 4.8.3"
