#ifndef CURVEFIT_H_
#define CURVEFIT_H_

#include <math.h>
#include <iostream>


using namespace std;
void fitExpFinetune(int n, const double *x, const double *y, double &a, double &b, double &c,bool tune_a,bool tune_b,bool tune_c);
void fitExp(int n, const double *x, const double *Y, double &a, double &b, double &c);
double predExp(double x,double a,double b,double c);
void fitExp_bfgs(int n, const double *x, const double *Y, double &a, double &b, double &c);

#endif
