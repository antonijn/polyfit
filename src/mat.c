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

#define EQ(x, y, e) (fabs(x - y) <= e)
#define NEQ(x, y, e) (!EQ(x, y, e))
#define EPS 0.0000001

/*
 * Elementary operations
 */

/* Swaps rows i and j in matrix m */
static void swaprows(struct mat m, int i, int j)
{
	for (int col = 0; col < m.numcols; ++col) {
		double temp = mget(m, i, col);
		mset(m, i, col, mget(m, j, col));
		mset(m, j, col, temp);
	}
}

/* Adds the scalar product of 'f' and row 'from' to row 'to' of matrix 'm' */
static void addrow(struct mat m, int to, int from, double f)
{
	assert(to != from || NEQ(f, -1.0, EPS));
	for (int j = 0; j < m.numcols; ++j)
		mset(m, to, j, mget(m, to, j) + f * mget(m, from, j));
}

/* Multiplies row i by f in matrix m */
static void mulrow(struct mat m, int i, double f)
{
	addrow(m, i, i, f - 1.0);
}

/*
 * Useful operations
 */

/* Inverts matrix 'w' from column 'fromc' */
static int invertfrom(struct mat w, int fromc)
{
	/* first we must get a one in w[fromc][fromc] */
	int rowtoswap = fromc;
	while (1) {
		if (NEQ(mget(w, rowtoswap, fromc), 0.0, EPS))
			break;
		if (rowtoswap >= w.numrows)
			return -1;
		++rowtoswap;
	}

	if (rowtoswap != fromc)
		swaprows(w, rowtoswap, fromc);

	mulrow(w, fromc, 1.0 / mget(w, fromc, fromc));

	/* now we elliminate */
	for (int i = 0; i < w.numrows; ++i) {
		if (i == fromc)
			continue;
		addrow(w, i, fromc, -mget(w, i, fromc));
	}

	return (fromc == w.numrows - 1) ? 0 : invertfrom(w, fromc + 1);
}

/* Matrix inversion */
int invert(struct mat in)
{
	assert(in.numrows == in.numcols);

	struct mat new;
	new.numrows = in.numrows;
	new.numcols = in.numcols * 2;
	new.rows = malloc(sizeof(double) * in.numrows * in.numcols * 2);
	for (int i = 0; i < in.numrows; ++i) {
		for (int j = 0; j < in.numcols; ++j) {
			mset(new, i, j, mget(in, i, j));
			mset(new, i, j + in.numcols, i == j);
		}
	}

	if (invertfrom(new, 0) < 0)
		return -1;

	for (int i = 0; i < in.numrows; ++i)
		for (int j = 0; j < in.numcols; ++j)
			mset(in, i, j, mget(new, i, j + in.numcols));

	return 0;
}

/*
 * Matrix multiplication
 */

/* Gets (l*r)[row][col] */
static double getmatmul(struct mat l, struct mat r, int row, int col)
{
	double res = 0.0;
	for (int k = 0; k < l.numcols; ++k)
		res += mget(l, row, k) * mget(r, k, col);
	return res;
}

/* Multiplies two matrices and stores result in the preallocated matrix 'to' */
void matmul(struct mat to, struct mat l, struct mat r)
{
	assert(l.numcols == r.numrows);
	assert(to.numrows == l.numrows);
	assert(to.numcols == r.numcols);

	for (int i = 0; i < to.numrows; ++i)
		for (int j = 0; j < to.numcols; ++j)
			mset(to, i, j, getmatmul(l, r, i, j));
}

/*
 * Matrix printing
 */

/* Prints row 'row' */
static void printrow(struct mat m, int row)
{
	assert(row >= 0);
	assert(row < m.numrows);

	printf("(");
	for (int i = 0; i < m.numcols; ++i) {
		double d = mget(m, row, i);
		if (d >= 0.0)
			putchar(' ');
		printf("%lf", d);
		printf((i == m.numcols - 1) ? ")" : ", ");
	}
}

/* Prints a matrix followed by a newline */
void printmat(struct mat m)
{
	for (int i = 0; i < m.numrows; ++i) {
		printrow(m, i);
		printf("\n");
	}
}

/* Prints a matrix multiplication */
void printmul(struct mat l, struct mat r, struct mat res)
{
	assert(l.numrows == r.numrows && r.numrows == res.numrows);

	for (int i = 0; i < l.numrows; ++i) {
		printrow(l, i);
		printrow(r, i);
		printf(i == l.numrows / 2 ? " = " : "   ");
		printrow(res, i);
		putchar('\n');
	}
}

