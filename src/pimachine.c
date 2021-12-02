#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define DIGITS_PER_ITERATION 14.1816474627254776555
#define BITS_PER_DIGIT 3.32192809488736234789 // log2(10)

typedef struct {
	mpz_t n;
	mpz_t d;
} mpz_frac;

void addFracArray(mpz_frac *rop, const mpz_frac arr[], unsigned long len) {
	mpz_t multiple, num;
	mpz_inits(multiple, num, NULL);

	// Store the first numerator of arr into rop to use as initial value
	mpz_set(rop->d, arr[0].d);

	// Find LCM of fraction array and store it as the denomenator
	for (unsigned long i = 1; i < len; i++) {
		mpz_lcm(rop->d, rop->d, arr[i].d);
	}

	/* Scale the numerator of each fraction to the new denomenator
	 * and sum it into the new numerator */
	for (unsigned long i = 0; i < len; i++) {
		mpz_div(multiple, rop->d, arr[i].d);
		mpz_mul(num, multiple, arr[i].n);
		mpz_add(rop->n, rop->n, num);
	}

	mpz_clears(multiple, num, NULL);
}

void calcSeries(mpz_frac *rop, unsigned long n) {
	mpz_t mnum, mden, l, x, k, m1, m2;
	mpz_frac facts[n];

	// Initialize the mpz_frac array
	for (unsigned long i = 0; i < n; i++)
		mpz_inits(facts[i].n, facts[i].d, NULL);

	// Initialize variables used in sum
	mpz_inits(m1, m2, NULL);

	// Initialize sequence variables with zero values
	mpz_init_set_si(mnum, 1);
	mpz_init_set_si(mden, 1);
	mpz_init_set_si(l, 13591409);
	mpz_init_set_si(x, 1);
	mpz_init_set_si(k, -6);

	// Initialize iteration and rop variable
	mpz_set(facts[0].n, l);
	mpz_set(facts[0].d, x);

	for (unsigned long q = 1; q < n; q++) {
		// Calculate k
		mpz_add_ui(k, k, 12);
		// Calculate m
		mpz_pow_ui(m1, k, 3);
		mpz_mul_ui(m2, k, 16);
		mpz_sub(m1, m1, m2);
		mpz_ui_pow_ui(m2, q, 3);
		mpz_mul(mnum, mnum, m1);
		mpz_mul(mden, mden, m2);
		// Calculate l
		mpz_add_ui(l, l, 545140134);
		// Calculate x
		mpz_mul_si(x, x, -262537412640768000L);

		// Derive sequence from m, l, x
		mpz_mul(facts[q].n, l, mnum);
		mpz_mul(facts[q].d, x, mden);
	}

	// Sum all of the generated fractions into rop.
	if (n > 1) {
		addFracArray(rop, facts, n);
	} else {
		mpz_set(rop->n, facts[0].n);
		mpz_set(rop->d, facts[0].d);
	}

	// Free the mpz_frac array
	for (unsigned long i = 0; i < n; i++)
		mpz_clears(facts[i].n, facts[i].d, NULL);

	// Free all of our variables
	mpz_clears(mnum, mden, l, x, k, m1, m2, NULL);
}

// WARNING: this returns an malloc'd string! You should probably free it later.
char *calcPi(unsigned long digits) {
	// Set precision of our float.
	unsigned long precisionBits = digits * BITS_PER_DIGIT + 1;
	unsigned long iterations = digits / DIGITS_PER_ITERATION + 1;
	mpf_set_default_prec(precisionBits);

	mpf_t c, pi, sumnf, sumdf, invSum;
	mpz_frac sum;
	mp_exp_t exp;
	mpf_inits(c, pi, sumnf, sumdf, invSum, NULL);
	mpz_inits(sum.n, sum.d, NULL);

	// Calculate constant C
	mpf_sqrt_ui(c, 10005);
	mpf_mul_ui(c, c, 426880);

	// Solve for pi
	calcSeries(&sum, iterations);
	mpf_set_z(sumnf, sum.n);
	mpf_set_z(sumdf, sum.d);
	mpf_div(invSum, sumdf, sumnf);
	mpf_mul(pi, invSum, c);

	// Generate string and free variables
	char *output = mpf_get_str(NULL, &exp, 10, digits, pi);
	mpf_clears(c, pi, sumnf, sumdf, invSum, NULL);
	mpz_clears(sum.n, sum.d, NULL);
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
