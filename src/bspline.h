#ifndef BSPLINE_H_
#define BSPLINE_H_

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
class cubicBspline{
    public:
    cubicBspline(int ncoeffs,double* breakpts);
    void get_xi(double x, double* bsplineX);
    private:
    int nbreak,ncoeffs;
    gsl_bspline_workspace *bw;
    gsl_vector *breakpts;
    gsl_vector *B;
};
#endif