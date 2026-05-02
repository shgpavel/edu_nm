	.text
	.align 1;

/* Frame slots, relative to ar6 after the prologue. */
#define S_BASE       0
#define S_EXP        1
#define S_MOD        2
#define S_K          3
#define S_M          4
#define S_R1         5
#define S_INV        6
#define S_M_INV      7
#define S_R_MOD      8
#define S_R2_MOD     9
#define S_BASE_MOD   10
#define S_BASE_MONT  11
#define S_RES_MONT   12
#define S_E_ODD      13
#define S_R2         14
#define S_MASK       15
#define S_B          16
#define S_R          17
#define S_E_P2       18
#define FRAME_WORDS  32

	.macro MUL_LO out, a, b, tmp
	gr7 = \b;
	\tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\out = gr7;
	.endm

	/* MUL_LO with one of the middle mul steps VLIW-paired with a MOVE
	 * (mid_dst = mid_src). pack_dst must be distinct from gr7, \tmp, \a, \b. */
	.macro MUL_LO_PACK out, a, b, tmp, mid_dst, mid_src
	gr7 = \b;
	\tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\mid_dst = \mid_src with \tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\out = gr7;
	.endm

	/* MUL_LO assuming gr7 already preloaded with \b. */
	.macro MUL_LO_PRELOADED out, a, tmp
	\tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\out = gr7;
	.endm

	/* MUL_LO assuming gr7 already preloaded with \b. Pairs one mid mul
	 * step with an inline MOVE. */
	.macro MUL_LO_PACK_PRELOADED out, a, tmp, mid_dst, mid_src
	\tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\mid_dst = \mid_src with \tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\out = gr7;
	.endm

	.macro UMUL_HILO a, b, lo, hi, tmp, ycopy, sign
	gr7 = \b;
	\ycopy = gr7 with \tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\lo = gr7 with \sign = \a A>> 31;
	\sign = \ycopy and \sign;
	\tmp = \tmp + \sign;
	\sign = \ycopy A>> 31;
	\hi = \a and \sign;
	\hi = \hi + \tmp;
	.endm

	/* UMUL_HILO with the final \hi = \hi + \tmp packed with a MOVE
	 * (tail_dst = tail_src). tail_dst must be distinct from \hi. */
	.macro UMUL_HILO_TAIL a, b, lo, hi, tmp, ycopy, sign, tail_dst, tail_src
	gr7 = \b;
	\ycopy = gr7 with \tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\lo = gr7 with \sign = \a A>> 31;
	\sign = \ycopy and \sign;
	\tmp = \tmp + \sign;
	\sign = \ycopy A>> 31;
	\hi = \a and \sign;
	\tail_dst = \tail_src with \hi = \hi + \tmp;
	.endm

	/* UMUL_HILO with one mid mul step paired with a MOVE
	 * (mid_dst = mid_src). mid_dst must be distinct from gr7, \tmp, \a, \b. */
	.macro UMUL_HILO_PACK a, b, lo, hi, tmp, ycopy, sign, mid_dst, mid_src
	gr7 = \b;
	\ycopy = gr7 with \tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\mid_dst = \mid_src with \tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\lo = gr7 with \sign = \a A>> 31;
	\sign = \ycopy and \sign;
	\tmp = \tmp + \sign;
	\sign = \ycopy A>> 31;
	\hi = \a and \sign;
	\hi = \hi + \tmp;
	.endm

	.macro UMUL_HILO_MNEG a, b, lo, hi, tmp, sign
	gr7 = \b;
	\sign = gr7 with \tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\lo = gr7 with \sign = \a A>> 31;
	\sign = \b and \sign;
	\tmp = \tmp + \sign;
	\hi = \a + \tmp;
	.endm

	.macro UMUL_HILO_SQ a, lo, hi, tmp, sign
	gr7 = \a;
	\tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\lo = gr7 with \sign = \a A>> 31;
	\sign = \a and \sign;
	\sign = \sign + \sign;
	\hi = \tmp + \sign;
	.endm

	/* UMUL_HILO_SQ with one mid mul step paired with a MOVE
	 * (mid_dst = mid_src). mid_dst must be distinct from gr7, \tmp, \a. */
	.macro UMUL_HILO_SQ_PACK a, lo, hi, tmp, sign, mid_dst, mid_src
	gr7 = \a;
	\tmp = \a *: gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\mid_dst = \mid_src with \tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\tmp = \a * gr7;
	\lo = gr7 with \sign = \a A>> 31;
	\sign = \a and \sign;
	\sign = \sign + \sign;
	\hi = \tmp + \sign;
	.endm

	.macro MONT_MUL a_slot, b_slot, out_slot, lbl
	gr0 = [ar6 + \a_slot];
	gr1 = [ar6 + \b_slot];
	UMUL_HILO gr0, gr1, gr2, gr3, gr4, gr5, gr6
	gr1 = [ar6 + S_M_INV];
	MUL_LO gr4, gr2, gr1, gr5
	gr1 = [ar6 + S_M];
	UMUL_HILO gr4, gr1, gr5, gr6, gr0, gr1, gr5
	gr4 = 0;
	gr4 = gr4 - gr2;
	gr4 = gr4 or gr2;
	gr4 = gr4 A>> 31;
	gr5 = 1;
	gr4 = gr4 and gr5;
	gr5 = 0;
	gr0 = gr3 + gr6;
	gr5 = gr5 + carry;
	gr0 = gr0 + gr4;
	gr5 = gr5 + carry;
	gr4 = 0;
	gr4 = gr4 - gr5;
	gr6 = [ar6 + S_R_MOD];
	gr6 = gr6 and gr4;
	gr0 = gr0 + gr6;
	gr5 = [ar6 + S_M];
	gr0 - gr5;
	if u< goto \lbl\()_store;
	gr0 = gr0 - gr5;
\lbl\()_store:
	[ar6 + \out_slot] = gr0;
	.endm

	.macro MONT_MUL_AR_PRE out_ar, lbl
	UMUL_HILO gr0, gr1, gr2, gr3, gr4, gr5, gr6
	\out_ar = gr2;
	MUL_LO_PACK gr4, gr2, ar4, gr5, gr1, ar3
	UMUL_HILO_TAIL gr4, gr1, gr2, gr6, gr0, gr1, gr5, gr4, \out_ar
	gr4 = gr4 + gr2;
	gr5 = 0 with gr4 = gr4 + carry;
	gr0 = gr3 + gr6;
	gr5 = gr5 + carry;
	gr0 = gr0 + gr4;
	gr5 = gr5 + carry;
	gr6 = ar5 with gr4 = -gr5;
	gr6 = gr6 and gr4;
	gr5 = ar3 with gr0 = gr0 + gr6;
	gr0 - gr5;
	if u< delayed goto \lbl\()_done;
	\out_ar = gr0;
	gr0 = gr0 - gr5;
	\out_ar = gr0;
\lbl\()_done:
	.endm

	.macro MONT_MUL_AR_PRE_MNEG out_ar, lbl
	UMUL_HILO gr0, gr1, gr2, gr3, gr4, gr5, gr6
	\out_ar = gr2;
	gr1 = ar4;
	MUL_LO_PACK gr4, gr2, gr1, gr5, gr1, ar3
	UMUL_HILO_MNEG gr4, gr1, gr2, gr6, gr0, gr5
	gr4 = \out_ar;
	gr4 = gr4 + gr2;
	gr5 = 0 with gr4 = gr4 + carry;
	gr0 = gr3 + gr6;
	gr5 = gr5 + carry;
	gr0 = gr0 + gr4;
	gr5 = gr5 + carry;
	gr6 = ar5 with gr4 = -gr5;
	gr6 = gr6 and gr4;
	gr5 = ar3 with gr0 = gr0 + gr6;
	gr0 - gr5;
	if u< delayed goto \lbl\()_done;
	\out_ar = gr0;
	gr0 = gr0 - gr5;
	\out_ar = gr0;
\lbl\()_done:
	.endm

	.macro MONT_SQUARE_AR_PRE out_ar, lbl
	UMUL_HILO_SQ_PACK gr0, gr2, gr3, gr4, gr5, gr1, ar4
	\out_ar = gr2 with gr7 = gr1;
	MUL_LO_PACK_PRELOADED gr4, gr2, gr5, gr1, ar3
	UMUL_HILO_TAIL gr4, gr1, gr2, gr6, gr0, gr1, gr5, gr4, \out_ar
	gr4 = gr4 + gr2;
	gr5 = 0 with gr4 = gr4 + carry;
	gr0 = gr3 + gr6;
	gr5 = gr5 + carry;
	gr0 = gr0 + gr4;
	gr5 = gr5 + carry;
	gr6 = ar5 with gr4 = -gr5;
	gr6 = gr6 and gr4;
	gr5 = ar3 with gr0 = gr0 + gr6;
	gr0 - gr5;
	if u< delayed goto \lbl\()_done;
	\out_ar = gr0;
	gr0 = gr0 - gr5;
	\out_ar = gr0;
\lbl\()_done:
	.endm

	.macro MONT_SQUARE_AR_PRE_MNEG out_ar, lbl
	UMUL_HILO_SQ gr0, gr2, gr3, gr4, gr5
	\out_ar = gr2;
	MUL_LO_PACK gr4, gr2, gr1, gr5, gr1, ar3
	UMUL_HILO_MNEG gr4, gr1, gr2, gr6, gr0, gr5
	gr4 = \out_ar;
	gr4 = gr4 + gr2;
	gr5 = 0 with gr4 = gr4 + carry;
	gr0 = gr3 + gr6;
	gr5 = gr5 + carry;
	gr0 = gr0 + gr4;
	gr5 = gr5 + carry;
	gr6 = ar5 with gr4 = -gr5;
	gr6 = gr6 and gr4;
	gr5 = ar3 with gr0 = gr0 + gr6;
	gr0 - gr5;
	if u< delayed goto \lbl\()_done;
	\out_ar = gr0;
	gr0 = gr0 - gr5;
	\out_ar = gr0;
\lbl\()_done:
	.endm

	.macro MONT_MUL_AR a_ar, b_ar, out_ar, lbl
	gr0 = \a_ar;
	gr1 = \b_ar;
	MONT_MUL_AR_PRE \out_ar, \lbl
	.endm

	.macro MONT_MUL_AR_MNEG a_ar, b_ar, out_ar, lbl
	gr0 = \a_ar;
	gr1 = \b_ar;
	MONT_MUL_AR_PRE_MNEG \out_ar, \lbl
	.endm

	/* MONT_REDUCE_AR: leaves the reduced word in gr6. Caller decides
	 * whether/where to store it. */
	.macro MONT_REDUCE_AR in_ar, lbl
	gr2 = \in_ar;
	gr1 = ar4;
	MUL_LO_PACK gr4, gr2, gr1, gr5, gr1, ar3
	UMUL_HILO gr4, gr1, gr5, gr6, gr0, gr1, gr3
	gr0 = 0;
	gr2;
	if =0 goto \lbl\()_no_carry;
	gr0 = 1;
\lbl\()_no_carry:
	gr6 = gr6 + gr0;
	.endm

	.globl _contest_pow_mod
	.type	_contest_pow_mod, @function
_contest_pow_mod:
	push ar6,gr6;
	push ar5,gr5;
	push ar4,gr4;
	push ar3,gr3;
	push ar2,gr2;
	push ar1,gr1;
	push ar0,gr0;
	ar6 = ar7;
	ar7 = ar7 + FRAME_WORDS;

	gr0 = [ar6 + -17];
	gr1 = [ar6 + -18];
	gr2 = [ar6 + -19];
	[ar6 + S_BASE] = gr0;             /* still needed by the pow2 path */
	[ar6 + S_EXP] = gr1;              /* still needed by the pow2 path */
	ar2 = gr0;                        /* stash base; ar2 is later reused for exp */

	gr3 = 1;
	gr2 - gr3;
	if <>0 delayed goto L_mod_not_one;
	gr1;
	vnul;
	gr7 = 0;
	goto L_return;

L_mod_not_one:
	if <>0 delayed goto L_exp_not_zero;
	gr0 = false;
	vnul;
	gr7 = 1;
	goto L_return;

L_exp_not_zero:
	/* gr2 holds modulus; gr0 = 0 (delay-slot clear); ar2 = base. */
	[ar6 + S_K] = gr0;

	/* gr2 already holds the modulus from the prologue (S_MOD = m initially). */
L_split_loop:
	gr1 = 1;
	gr0 = gr2 and gr1;
	gr0;
	if <>0 delayed goto L_split_done_loaded;
	vnul;
	gr3 = [ar6 + S_K];
	gr3 = gr3 + 1;
	[ar6 + S_K] = gr3;
	gr2 = gr2 >> 1;
	goto L_split_loop;

L_split_done_loaded:
	/* gr2 = m (odd part). m is now ar3-resident; no consumer reads S_M. */
	ar3 = gr2;                        /* keep m register-resident */
	gr1 = 1;
	gr1 - gr2;
	if u>= goto L_pow2_start;

	/* inv = m^(-1) mod 2^32 by four Newton lifts. m stays in gr0. */
	gr0 = gr2;
	gr1 = gr0;
	MUL_LO gr2, gr0, gr1, gr3
	gr4 = 2;
	gr2 = gr4 - gr2;
	MUL_LO gr1, gr1, gr2, gr3
	MUL_LO gr2, gr0, gr1, gr3
	gr4 = 2;
	gr2 = gr4 - gr2;
	MUL_LO gr1, gr1, gr2, gr3
	MUL_LO gr2, gr0, gr1, gr3
	gr4 = 2;
	gr2 = gr4 - gr2;
	MUL_LO gr1, gr1, gr2, gr3
	MUL_LO gr2, gr0, gr1, gr3
	gr4 = 2;
	gr2 = gr4 - gr2;
	MUL_LO gr1, gr1, gr2, gr3
	gr2 = false;
	gr2 = gr2 - gr1;
	ar4 = gr2;                        /* m_inv register-resident */

	/* r_mod = 2^32 mod m = (-m) % m. gr0 still holds m. */
	gr1 = false;
	gr1 = gr1 - gr0;
	gr1 - gr0;
	if u< delayed goto L_r_mod_inline;
	vnul;
	vnul;
	delayed call UMod32;
	push gr1;
	push gr0;
	goto L_r_mod_set;
L_r_mod_inline:
	gr7 = gr1;                        /* (-m) when (-m) < m */
L_r_mod_set:
	ar0 = gr7;                        /* Montgomery representation of 1 = R_MOD */
	ar5 = gr7;

	/* base_mod = base % m. */
	gr0 = ar2;                        /* base */
	gr1 = ar3;                        /* m */
	gr0 - gr1;
	if u< delayed goto L_base_mod_done;
	gr2 = gr0 - gr1;                  /* gr2 = candidate base_mod */
	vnul;
	gr2 - gr1;
	if u< delayed goto L_base_mod_have_gr2;
	vnul;
	vnul;
	delayed call UMod32;
	push gr0;
	push gr1;
	gr2 = gr7;
L_base_mod_have_gr2:
	gr0 = gr2;                        /* base_mod -> gr0 for the next step */
L_base_mod_done:
	/* gr0 = base_mod. Build base_mont = (base_mod << 32) % m via UMod64. */
	gr1 = gr0;                        /* low half of dividend = 0, high = base_mod */
	gr0 = false;
	gr2 = ar3;                        /* m */
	gr3 = false;
	push gr0;
	push gr1;
	push gr2;
	push gr3;
	delayed call UMod64;
	vnul;
	vnul;
	pop gr2;
	pop gr1;
	ar7 -= 2;
	ar1 = gr1;                        /* base_mont */
	gr0 = [ar6 + S_EXP];
	ar2 = gr0;
	/* ar0 = R_MOD already set above; this is Montgomery representation of 1. */

L_odd_exp_loop:
	gr0 = ar2;
	gr0;
L_odd_exp_loop_tested:
	if =0 delayed goto L_odd_exp_done;
	gr1 = 1;
	gr2 = gr0 and gr1;
	gr2;
	if =0 delayed goto L_odd_skip_mul;
	gr0 = gr0 >> 1;
	ar2 = gr0;
	MONT_MUL_AR ar0, ar1, ar0, L_mont_res
L_odd_skip_mul:
	gr0 = ar2;
	gr0;
	if =0 delayed goto L_odd_exp_done;
	gr0 = ar1;
	MONT_SQUARE_AR_PRE ar1, L_mont_square
	delayed goto L_odd_exp_loop_tested;
	gr0 = ar2;
	gr0;

L_odd_exp_done:
	MONT_REDUCE_AR ar0, L_mont_final
	/* gr6 = R1 (the demontgomeryized residue). */
	gr0 = [ar6 + S_K];
	gr0;
	if =0 delayed goto L_return;
	gr7 = gr6;                        /* delay slot: R1 is the result */
	[ar6 + S_R1] = gr6;               /* k != 0: spill for CRT path */
	goto L_pow2_start;

L_pow2_start:
	/* mask = (1 << k) - 1, built without a variable shift. */
	gr0 = [ar6 + S_K];
	gr1 = 0;
	gr2 = gr0;
L_mask_loop:
	gr1 = gr1 + gr1;
	gr3 = 1;
	gr1 = gr1 or gr3;
	gr2 = gr2 - 1;
	gr2;
	if <>0 goto L_mask_loop;
	[ar6 + S_MASK] = gr1;

	gr0 = [ar6 + S_BASE];
	gr0 = gr0 and gr1;
	[ar6 + S_B] = gr0;
	gr0 = 1;
	[ar6 + S_R] = gr0;
	gr0 = [ar6 + S_EXP];
	[ar6 + S_E_P2] = gr0;

L_pow2_loop:
	gr0 = [ar6 + S_E_P2];
	gr0;
	if =0 goto L_pow2_store;
	gr1 = 1;
	gr2 = gr0 and gr1;
	gr2;
	if =0 goto L_pow2_skip_mul;
	gr0 = [ar6 + S_R];
	gr1 = [ar6 + S_B];
	MUL_LO gr2, gr0, gr1, gr3
	gr4 = [ar6 + S_MASK];
	gr2 = gr2 and gr4;
	[ar6 + S_R] = gr2;
L_pow2_skip_mul:
	gr0 = [ar6 + S_E_P2];
	gr0 = gr0 >> 1;
	[ar6 + S_E_P2] = gr0;
	gr0;
	if =0 goto L_pow2_loop;
	gr0 = [ar6 + S_B];
	MUL_LO gr2, gr0, gr0, gr3
	gr4 = [ar6 + S_MASK];
	gr2 = gr2 and gr4;
	[ar6 + S_B] = gr2;
	goto L_pow2_loop;

L_pow2_store:
	gr0 = [ar6 + S_R];
	[ar6 + S_R2] = gr0;

L_pow2_done:
	gr0 = ar3;                        /* m */
	gr1 = 1;
	gr0 - gr1;
	if <>0 goto L_crt;
	gr7 = [ar6 + S_R2];
	goto L_return;

L_crt:
	gr0 = [ar6 + S_R2];
	gr1 = [ar6 + S_R1];
	gr2 = gr0 - gr1;
	gr3 = [ar6 + S_MASK];
	gr2 = gr2 and gr3;
	gr0 = ar4;                        /* m_inv */
	gr5 = 0;
	gr0 = gr5 - gr0;
	gr0 = gr0 and gr3;
	MUL_LO gr4, gr2, gr0, gr5
	gr4 = gr4 and gr3;
	gr0 = ar3;                        /* m */
	MUL_LO gr2, gr0, gr4, gr5
	gr1 = [ar6 + S_R1];
	gr7 = gr1 + gr2;

L_return:
	ar7 = ar7 + -FRAME_WORDS;
	pop ar0,gr0;
	pop ar1,gr1;
	pop ar2,gr2;
	pop ar3,gr3;
	pop ar4,gr4;
	pop ar5,gr5;
	pop ar6,gr6;
	return;
	.size	_contest_pow_mod, .-_contest_pow_mod
