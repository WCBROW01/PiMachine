#include <cstdio>
#include <cstdlib>
#include <boost/multiprecision/gmp.hpp>

#define DIGITS_PER_ITERATION 14.1816474627254776555
#define BITS_PER_DIGIT 3.32192809488736234789 // log2(10)
#define FRACARRAY_MAXLEN 4

using namespace boost::multiprecision;

// Fractions will be initialized to 0/1 by default.
typedef struct {
	mpz_int n = 0;
	mpz_int d = 1;
} mpz_frac;

void addFracArray(mpz_frac &rop, const mpz_frac arr[], unsigned long len) {
	mpz_frac tempFrac;


	// Find LCM of fraction array and store it as the denomenator
	for (unsigned long i = 0; i < len; i++) {
		tempFrac.d = lcm(tempFrac.d, arr[i].d);
	}

	/* Scale the numerator of each fraction to the new denomenator
	 * and sum it into the new numerator */
	for (unsigned long i = 0; i < len; i++) {
		tempFrac.n += (tempFrac.d / arr[i].d) * arr[i].n;
	}

	// Merge fractions
	mpz_frac ropNew = rop;
	ropNew.d = lcm(tempFrac.d, rop.d);
	if (ropNew.n == 0) ropNew.n = tempFrac.n;
	else ropNew.n += (tempFrac.d / rop.d) * rop.n;
	rop = ropNew;
}

void calcSeries(mpz_frac &rop, unsigned long n) {
	mpz_frac facts[FRACARRAY_MAXLEN];

	// Initialize sequence variables with zero values
	mpz_int mnum = 1;
	mpz_int mden = 1;
	mpz_int l = 13591409;
	mpz_int x = 1;
	mpz_int k = -6;

	// Initialize iteration and rop variable
	facts[0].n = l;
	facts[0].d = x;

	for (unsigned long q = 1; q < n; q++) {
		// Calculate k
		k += 12;
		// Calculate m
		mnum *= pow(k, 3) - (k * 16);
		mden *= pow((mpz_int) q, 3);
		// Calculate l
		l += 545140134;
		// Calculate x
		x *= -262537412640768000L;

		// Derive sequence from m, l, x
		facts[q % FRACARRAY_MAXLEN].n = l * mnum;
		facts[q % FRACARRAY_MAXLEN].d = x * mden;

		if ((q + 1) % FRACARRAY_MAXLEN == 0) {
			addFracArray(rop, facts, FRACARRAY_MAXLEN);
		}
	}

	// Sum all of the generated fractions into rop.
	if (n % FRACARRAY_MAXLEN > 0)
		addFracArray(rop, facts, n % FRACARRAY_MAXLEN);
}

// WARNING: this returns an malloc'd string! You should probably free it later.
char *calcPi(unsigned long digits) {
	// Set precision of our float.
	unsigned long precisionBits = digits * BITS_PER_DIGIT + 1;
	unsigned long iterations = digits / DIGITS_PER_ITERATION + 1;
	mpf_float::default_precision(precisionBits);

	mpf_float pi, invSum;
	mpz_frac sum;
	mp_exp_t exp;

	// Calculate constant C
	mpf_float c = sqrt(10005);
	c *= 426880;

	// Solve for pi
	calcSeries(sum, iterations);
	invSum = (mpf_float) sum.d / sum.n;
	pi = invSum * c;

	// Generate string and free variables
	char *output = mpf_get_str(NULL, &exp, 10, digits, pi.backend().data());
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
