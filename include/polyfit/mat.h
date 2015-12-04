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

#ifndef MAT_H
#define MAT_H

#include <stdlib.h>

struct point {
	double x, y;
};

struct mat {
	double *rows;
	size_t numrows, numcols;
};

void printmat(struct mat m);
void printmul(struct mat l, struct mat r, struct mat res);

int invert(struct mat m);
void matmul(struct mat to, struct mat l, struct mat r);

static inline double mget(struct mat m, int row, int col)
{
	return m.rows[col + row*m.numcols];
}

static inline double mset(struct mat m, int row, int col, double v)
{
	return (m.rows[col + row*m.numcols] = v);
}

#endif

