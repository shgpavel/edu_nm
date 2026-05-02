	.text
	.align 1
	.globl _contest_a_plus_b
	.type _contest_a_plus_b, @function
_contest_a_plus_b:
	gr7 = [ar7 + -3];
	gr0 = [ar7 + -4];
	gr7 = gr7 + gr0;
	return;
	.size _contest_a_plus_b, .-_contest_a_plus_b
