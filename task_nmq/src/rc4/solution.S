	.text
	.align 1
	.globl _contest_rc4_init
	.type _contest_rc4_init, @function
_contest_rc4_init:
	push ar6,gr6;
	push ar5,gr5;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	ar6 = ar7;

	ar0 = ar6 + -16;
	gr2 = [ar0 + -1];	// ctx
	gr4 = [ar0 + -2];	// key
	gr5 = [ar0 + -3];	// key_len

	// Args are fully read; reuse ar6 as the ctx base.
	ar6 = gr2;

	// --- Fill S[i] = i for i = 0..255 ---
	ar0 = ar6;		// advancing S ptr
	gr0 = 0;
	gr1 = 255;
	gr1;
L_fill:
	if <>0 delayed goto L_fill;
	[ar0++] = gr0 with gr0 = gr0 + 1;
	gr1 = gr1 - 1;

	// --- KSA ---
	// ar6 = ctx base (permanent; j in gr6 makes [ar6 + gr6] valid)
	// ar2 = advancing key ptr, gr4 = key start, gr5 = key_len
	// ar0 = &S[i] advancing ptr, gr1 = countdown, gr3 = 255
	ar2 = gr4;		// key ptr
	ar0 = ar6;		// &S[0]
	gr6 = 0;		// j
	gr3 = 255;
	gr1 = 256;
	gr0 = gr5;		// key countdown
L_ksa:
	gr7 = [ar0];
	gr2 = [ar2] with gr6 = gr7 + gr6;
	ar2 = ar2 + 1 with gr6 = gr2 + gr6;
	gr6 = gr6 and gr3;
	gr0 = gr0 - 1;
	if =0 delayed goto L_ksa_key_reset;
	gr2 = [ar6 + gr6];
	gr1 = gr1 - 1;
	if <>0 delayed goto L_ksa;
	[ar6 + gr6] = gr7;
	[ar0++] = gr2;
	goto L_ksa_done;
L_ksa_key_reset:
	ar2 = gr4;
	if <>0 delayed goto L_ksa;
	[ar6 + gr6] = gr7;
	[ar0++] = gr2 with gr0 = gr5;
L_ksa_done:

	// Store i=0, j=0 into ctx->i and ctx->j.
	ar3 = ar6 + 256;
	gr0 = 0;
	[ar3] = gr0;
	ar2 = ar6 + 257;
	[ar2] = gr0;

	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size _contest_rc4_init, .-_contest_rc4_init

	.align 1
	.globl _contest_rc4_crypt
	.type _contest_rc4_crypt, @function
_contest_rc4_crypt:
	push ar6,gr6;
	push ar5,gr5;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;

	ar7 = ar7 + 8;
	ar6 = ar7 + -8;
	ar0 = ar6 + -16;

	gr0 = [ar0 + -1];	// ctx
	gr2 = [ar0 + -2];	// input
	ar1 = [ar0 + -3];	// output (permanent)
	gr4 = [ar0 + -4];	// len

	ar4 = gr2;		// input ptr (permanent)

	// Read i and j from ctx
	ar0 = gr0;
	ar3 = ar0 + 256;	// &ctx->i
	gr2 = [ar3];		// i

	ar0 = ar0 + 257;	// &ctx->j
	gr5 = [ar0];		// j (kept in gr5 throughout loop)

	// Set up 3 permanent ctx copies for ar+gr addressing:
	//   ar2=ctx, gr2=i -> [ar2 + gr2]
	//   ar5=ctx, gr5=j -> [ar5 + gr5]
	//   ar6=ctx, gr6=t -> [ar6 + gr6]
	// ar6 frame ptr no longer needed after args are read.
	ar2 = gr0;
	ar5 = gr0;
	ar6 = gr0;
	gr3 = 255;

	gr4;
	if v<= delayed goto L_crypt_done;
	gr2 = gr2 + 1;
	gr2 = gr2 and gr3;

L_crypt_loop:
	gr7 = [ar2 + gr2];
	gr1 = [ar4++] with gr5 = gr7 + gr5;
	gr5 = gr5 and gr3;
	gr6 = [ar5 + gr5] with gr4 = gr4 - 1;
	nul;
	[ar5 + gr5] = gr7 with gr7 = gr7 + gr6 noflags;
	[ar2 + gr2] = gr6 with gr6 = gr7 and gr3 noflags;
	gr7 = [ar6 + gr6] with gr2 = gr2 + 1 noflags;
	if <>0 delayed goto L_crypt_loop;
	gr0 = gr7 xor gr1;
	[ar1++] = gr0 with gr2 = gr2 and gr3;

L_crypt_done:
	// Write i,j back via ar2=ctx (frame ptr ar6 was reused; ar2 is permanent)
	ar0 = ar2 + 256;
	[ar0] = gr2;
	ar0 = ar2 + 257;
	[ar0] = gr5;

	ar7 = ar7 + -8;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size _contest_rc4_crypt, .-_contest_rc4_crypt
