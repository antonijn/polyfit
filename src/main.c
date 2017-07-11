/* Copyright (c) 2015 antonijn
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <polyfit/mat.h>

/* Calculates x^n */
static double npow(double x, int n)
{
	assert(n >= 0);
	if (n == 0)
		return 1;

	return x * npow(x, n - 1);
}

/* Nicely prints "{coef} x^{exp}" */
static void printxexp(double coef, int exp)
{
	printf("%lf", coef);
	switch (exp) {
	case 0:
		break;
	case 1:
		printf(" x");
		break;
	default:
		printf(" x^%d", exp);
	}
}

/* Prints the values of column vector 'm' as coefficients of a polynomial */
static void polyprint(struct mat m)
{
	if (m.numrows == 0)
		return;

	printf("P(x) = ");

	int fst = -1;
	while (fabs(m.rows[++fst]) < 0.000001)
		;

	for (int i = fst; i < m.numrows; ++i) {
		if (fabs(m.rows[i]) < 0.000001)
			continue;

		if (i != fst || m.rows[i] < 0)
			printf(" %c ", (m.rows[i] >= 0) ? '+' : '-');
		printxexp(fabs(m.rows[i]), (int)m.numrows - (i + 1));
	}
}

int main(int argc, char **argv)
{
	if (argc == 1)
		return 0;

	const int matsize = argc - 1;
	struct mat xs;
	xs.numcols = matsize;
	xs.numrows = matsize;
	xs.rows = malloc(sizeof(double) * matsize * matsize);

	struct mat ys;
	ys.numcols = 1;
	ys.numrows = matsize;
	ys.rows = malloc(sizeof(double) * matsize);

	for (int i = 0; i < matsize; ++i) {
		char *pt = argv[i + 1];
		double x, y;
		if (sscanf(pt, "%lf,%lf", &x, &y) != 2) {
			fprintf(stderr, "Invalid input\n");
			return EXIT_FAILURE;
		}

		for (int n = 0; n < matsize; ++n)
			mset(xs, i, matsize - (n + 1), npow(x, n));

		mset(ys, i, 0, y);
	}

	if (invert(xs) < 0)
		return EXIT_FAILURE;

	struct mat to;
	to.numcols = 1;
	to.numrows = matsize;
	to.rows = malloc(sizeof(double) * matsize);
	matmul(to, xs, ys);

	polyprint(to);
	printf("\n");

	free(xs.rows);
	free(ys.rows);
	free(to.rows);
	return 0;
}

