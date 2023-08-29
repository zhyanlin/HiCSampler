/*
 * optimize.h
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

#include "utils.h"

double greedyOptimizeT(const double &o, double t, const double allNei[], const int &nNei, double bias);
double evaluateFullLikelihoodPatch(int ii, int jj);
double evaluateFullLikelihoodPatch(int ii, int jj, double newVal);
double evaluateFullLikelihood();
double evaluateFullLikelihood2();
double evaluateLikelihood(const float &currentT, const int &obs, const float &medNei, const float &bias);
int getNeighbors(const int &i, const int &j, int &nNei, double allNei[]);

#endif /* OPTIMIZE_H_ */
