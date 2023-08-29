#include "bspline.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
using namespace std;

/* number of data points to fit */
#define N 200

/* number of fit coefficients */

/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
#define NBREAK (NCOEFFS - 2)

cubicBspline::cubicBspline(int ncoeffs , double *breakpts)
{
    int k = 4;
    nbreak = ncoeffs + 2 - k;
    this->ncoeffs = ncoeffs;
    bw = gsl_bspline_alloc(k, nbreak);
    this->breakpts = gsl_vector_alloc(nbreak);
    for (int i = 0; i < nbreak; ++i)
        gsl_vector_set(this->breakpts, i, breakpts[i]);
    gsl_bspline_knots(this->breakpts, bw);
    B = gsl_vector_alloc(ncoeffs);
};
int cubicBspline::get_xi(double xi, double *bsplineX)
{
    gsl_bspline_eval(xi, B, bw);
    for (int j = 0; j < ncoeffs; ++j)
        bsplineX[j] = gsl_vector_get(B, j);
};
