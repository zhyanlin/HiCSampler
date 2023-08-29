/*
 * utils.h
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <math.h>
#include <zlib.h>
#include <pthread.h>
#include "inlineFunc.h"
#include "contactMap.h"

#define SAFE 0

using namespace std;

#define MAX_NBNONMAX 10
#define BIAS_PSEUDOCOUNT 5
#define MAX_BIAS 10

#define NORM_RES 10000
#define NFACTORS 13

extern float factor[NFACTORS];
extern float logFactor[NFACTORS];
extern int fullRange;
extern float *logKFact;
extern float *bias, *piror_k, *piror_theta;
extern double *precomputedNormLogPDF;
extern int floatFullRange;
extern int chr1, chr2;

double norm_Randn(float mu, float sigma);

float normpdf(const float &x, const float &m, const float &s);

// precompute log k factorial
void precomputeLogKFact();
void precomputeStuff();

double evaluateAgainstTest();
double evaluateAgainstTest(UpperDiag<float> *currentT);
float getMedianNeighbor(int i, int j);
void setBoundaries();

void computeTMatrix_fixed();
void computeBias();
void readFullOMatrix(char *fn, UpperDiag<float> *mat, bool setBias);
void getMatrixSize(char *fn);
void allocateMemory();

float getMedianNeighbor(int i, int j);
void outputSparseMatrix(char *fn);
void outputSparseBoundaryMatrix(char *fn);
void outputSparseMatrixGZ(char *fn);
double median(int nNei, float *allNei);

struct annotation{
	double *len;
	double *gc;
	double *map;
};
extern annotation anno;
#endif /* UTILS_H_ */
