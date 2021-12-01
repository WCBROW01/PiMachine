#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define DIGITS_PER_ITERATION 14.1816474627254776555
#define BITS_PER_DIGIT 3.32192809488736234789 // log2(10)

void calcSeries(mpz_t ropn, mpz_t ropd, unsigned long n) {
	mpz_t numz, m, l, x, k, m1, m2;

	// Initialize variables used in sum
	mpz_inits(numz, m1, m2, NULL);

	// Initialize sequence variables with zero values
	mpz_init_set_si(m, 1);
	mpz_init_set_si(l, 13591409);
	mpz_init_set_si(x, 1);
	mpz_init_set_si(k, -6);

	// Derive numerator from m and l (first iteration)
	mpz_mul(numz, m, l);

	// Add this iteration of the sequence to the series
	mpz_add(ropn, ropn, numz);
	mpz_add(ropd, ropd, x);

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

		// Derive numerator from m, l
		mpz_mul(numz, m, l);

		// Add this iteration of the sequence to the series
		mpz_add(ropn, ropn, numz);
		mpz_add(ropd, ropd, x);
	}

	// Free all of our variables
	mpz_clears(numz, m, l, x, k, m1, m2, NULL);
}

// WARNING: this returns an malloc'd string! You should probably free it later.
char *calcPi(unsigned long digits) {
	// Set precision of our float.
	unsigned long precisionBits = (digits * BITS_PER_DIGIT) + 1;
	mpf_set_default_prec(precisionBits);

	mpf_t c, pi, invSum, sumnf, sumdf;
	mpz_t sumnz, sumdz;
	mp_exp_t exp;
	mpf_inits(c, pi, invSum, sumnf, sumdf, NULL);
	mpz_inits(sumnz, sumdz, NULL);

	// Calculate constant C
	mpf_sqrt_ui(c, 10005);
	mpf_mul_ui(c, c, 426880);

	// Solve for the series
	calcSeries(sumnz, sumdz, digits / DIGITS_PER_ITERATION + 1);

	// Convert the numerator and denomenator to floats
	mpf_set_z(sumnf, sumnz);
	mpf_set_z(sumdf, sumdz);

	// Get the inverse of the sum
	mpf_div(invSum, sumdf, sumnf);

	// Multiply by c to get pi
	mpf_mul(pi, invSum, c);

	// Generate string and free variables
	char *output = mpf_get_str(NULL, &exp, 10, digits, pi);
	printf("%lu\n", exp);
	mpf_clears(c, pi, invSum, sumnf, sumdf, NULL);
	mpz_clears(sumnz, sumdz, NULL);
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
