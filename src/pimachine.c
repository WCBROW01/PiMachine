#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define DIGITS_PER_ITERATION 14.1816474627254776555
#define BITS_PER_DIGIT 3.32192809488736234789 // log2(10)

void calcSeries(mpf_t rop, unsigned long n) {
	mpz_t numz, m, l, x, k, m1, m2;
	mpf_t numf, denf, iteration;

	// Initialize variables used in sum
	mpz_inits(numz, m1, m2, NULL);
	mpf_inits(numf, denf, iteration, NULL);

	// Initialize sequence variables with zero values
	mpz_init_set_si(m, 1);
	mpz_init_set_si(l, 13591409);
	mpz_init_set_si(x, 1);
	mpz_init_set_si(k, -6);

	// Derive sequence from m, l, x (first iteration)
	mpz_mul(numz, m, l);
	mpf_set_z(numf, numz);
	mpf_set_z(denf, x);
	mpf_div(iteration, numf, denf);

	// Add this iteration of the sequence to the series
	mpf_add(rop, rop, iteration);

	for (unsigned long q = 1; q < n; q++) {
		// Calculate k
		mpz_add_ui(k, k, 12);
		// Calculate m
		mpz_pow_ui(m1, k, 3);
		mpz_mul_ui(m2, k, 16);
		mpz_sub(m1, m1, m2);
		mpz_cdiv_q_ui(m1, m1, q * q * q);
		mpz_mul(m, m, m1);
		// Calculate l
		mpz_add_ui(l, l, 545140134);
		// Calculate x
		mpz_mul_si(x, x, -262537412640768000L);

		// Derive sequence from m, l, x
		mpz_mul(numz, m, l);
		mpf_set_z(numf, numz);
		mpf_set_z(denf, x);
		mpf_div(iteration, numf, denf);

		// Add this iteration of the sequence to the series
		mpf_add(rop, rop, iteration);
	}

	// Free all of our variables
	mpz_clears(numz, m, l, x, k, m1, m2, NULL);
	mpf_clears(numf, denf, iteration, NULL);
}

// WARNING: this returns an malloc'd string! You should probably free it later.
char *calcPi(unsigned long digits) {
	// Set precision of our float.
	unsigned long precisionBits = (digits * BITS_PER_DIGIT) + 1;
	mpf_set_default_prec(precisionBits);

	mpf_t c, pi, sum;
	mp_exp_t exp;
	mpf_inits(c, pi, sum, NULL);

	// Calculate constant C
	mpf_sqrt_ui(c, 10005);
	mpf_mul_ui(c, c, 426880);

	// Solve for pi
	calcSeries(sum, digits / DIGITS_PER_ITERATION + 1);
	mpf_ui_div(sum, 1, sum);
	mpf_mul(pi, sum, c);

	// Generate string and free variables
	char *output = mpf_get_str(NULL, &exp, 10, digits, pi);
	mpf_clears(c, pi, sum, NULL);
	return output;
}

int main(int argc, char *argv[]) {
	// Will add more arg possibilities.
	switch(argc) {
		case 1:
			printf("Please provide the number of digits you want to calculate.");
			return 1;
	}
	char *piStr = calcPi(strtoul(argv[1], NULL, 10));
	printf("%.1s.%s\n", piStr, piStr + 1);
	free(piStr);
	return 0;
}
