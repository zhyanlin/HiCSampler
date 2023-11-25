#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int n_threads;
extern int option_maxFragDistance;
extern int k;
extern bool option_diffusion;
extern int option_mcmcMaxIter;
extern bool method_mrf;
extern bool method_fixed;
extern bool method_kde;
extern bool option_lognormal;
extern bool option_median;
extern double option_kdeBandwidth;
extern double option_kdeMaxBandwidth;
extern double option_minSigma;
extern double option_sigmaMultiplier;
extern bool option_kdeInitialization;
extern int option_diagBand;
extern int option_kdeMinCount;
extern double option_boundaryDensity;
extern double option_boundaryDenisty;
extern double option_sigma,option_w;
extern double option_averageFragmentSize;
extern bool option_inputFullMatrix;
extern bool option_inputSparseMatrix;
extern bool option_noBias;
extern double option_pseudocount;
extern bool option_mcmc;
extern bool option_gibbs;
extern double option_boundaryKS;
extern bool option_diagRenormalization;
extern bool option_outputNormalized;
extern double option_minOutput;
extern int option_mrfMaxIter;
extern int option_firstRow;
extern int option_firstCol;
extern int option_lastRow;
extern int option_lastCol;
extern int Tr, Tc;
extern char  option_prior[100];
extern bool option_outputSymmetric;
extern char option_initializationMatrix[1000];
extern char option_testMatrix[1000];
extern char option_bias[1000];
extern char option_prior_file[1000];
extern char option_annotationFile[1000];
extern char option_boundaryOutput[1000];
extern int bandSize; // number of fragments to consider in terms of maxi genomic distance when run it under a band configuration
extern int stepSize;
extern bool outputSamples;
void readOptions(int argc, char *argv[]);
