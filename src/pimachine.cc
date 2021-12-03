#include <cstdio>
#include <cstdlib>
#include <boost/multiprecision/gmp.hpp>

#define DIGITS_PER_ITERATION 14.1816474627254776555
#define BITS_PER_DIGIT 3.32192809488736234789 // log2(10)

using namespace boost::multiprecision;

typedef struct {
	mpz_int n;
	mpz_int d;
} mpz_frac;

void addFracArray(mpz_frac *rop, const mpz_frac arr[], unsigned long len) {
	mpz_int multiple, num;
	//mpz_inits(multiple, num, NULL);

	// Store the first numerator of arr into rop to use as initial value
	rop->d = arr[0].d;

	// Find LCM of fraction array and store it as the denomenator
	for (unsigned long i = 1; i < len; i++) {
		rop->d = lcm(rop->d, arr[i].d);
	}

	/* Scale the numerator of each fraction to the new denomenator
	 * and sum it into the new numerator */
	for (unsigned long i = 0; i < len; i++) {
		multiple = rop->d / arr[i].d;
		num = multiple * arr[i].n;
		rop->n += num;
	}

	//mpz_clears(multiple, num, NULL);
}

void calcSeries(mpz_frac *rop, unsigned long n) {
	mpz_int mnum, mden, l, x, k, m1, m2;
	mpz_frac facts[n];

	// Initialize the mpz_frac array
	//for (unsigned long i = 0; i < n; i++)
		//mpz_inits(facts[i].n, facts[i].d, NULL);

	// Initialize variables used in sum
	//mpz_inits(m1, m2, NULL);

	// Initialize sequence variables with zero values
	mnum = 1;
	mden = 1;
	l = 13591409;
	x = 1;
	k = -6;

	// Initialize iteration and rop variable
	facts[0].n = l;
	facts[0].d = x;

	for (unsigned long q = 1; q < n; q++) {
		// Calculate k
		k += 12;
		// Calculate m
		m1 = pow(k, 3);
		m2 = k * 16;
		m1 -= m2;
		m2 = q * q * q;
		mnum *= m1;
		mden *= m2;
		// Calculate l
		l += 545140134;
		// Calculate x
		x *= -262537412640768000L;

		// Derive sequence from m, l, x
		facts[q].n = l * mnum;
		facts[q].d = x * mden;
	}

	// Sum all of the generated fractions into rop.
	if (n > 1) {
		addFracArray(rop, facts, n);
	} else {
		rop->n = facts[0].n;
		rop->d = facts[0].d;
	}

	// Free the mpz_frac array
	//for (unsigned long i = 0; i < n; i++)
		//mpz_clears(facts[i].n, facts[i].d, NULL);

	// Free all of our variables
	//mpz_clears(mnum, mden, l, x, k, m1, m2, NULL);
}

// WARNING: this returns an malloc'd string! You should probably free it later.
char *calcPi(unsigned long digits) {
	// Set precision of our float.
	unsigned long precisionBits = digits * BITS_PER_DIGIT + 1;
	unsigned long iterations = digits / DIGITS_PER_ITERATION + 1;
	mpf_float::default_precision(precisionBits);

	mpf_float c, pi, sumnf, sumdf, invSum;
	mpz_frac sum;
	mp_exp_t exp;
	//mpf_inits(c, pi, sumnf, sumdf, invSum, NULL);
	//mpz_inits(sum.n, sum.d, NULL);

	// Calculate constant C
	c = sqrt(10005);
	c *= 426880;

	// Solve for pi
	calcSeries(&sum, iterations);
	//mpf_set_z(sumnf, sum.n.backend().data());
	//mpf_set_z(sumdf, sum.d.backend().data());
	sumnf = (mpf_float) sum.n;
	sumdf = (mpf_float) sum.d;
	invSum = sumdf / sumnf;
	pi = invSum * c;

	// Generate string and free variables
	char *output = mpf_get_str(NULL, &exp, 10, digits, pi.backend().data());
	//mpf_clears(c, pi, sumnf, sumdf, invSum, NULL);
	//mpz_clears(sum.n, sum.d, NULL);
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
