/*
 * optimize.cpp
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#include "optimize.h"
#include "contactMap.h"
#include "readOptions.h"

int getNeighbors(const int &i, const int &j, int &nNei, double allNei[])
{
  nNei = 0;
  float sumNeighbor = 0;
  int ii, jj, off_j;
  int startRow = O.startRow(), endRow = O.endRow();
  for (int off_i = -1; off_i <= 1; off_i++)
  {
    for (off_j = -1; off_j <= 1; off_j++)
    {
      if (off_i == 0 && off_j == 0)
        continue;

      ii = i + off_i;
      jj = j + off_j;
      if (ii < jj && ii >= startRow && ii < endRow && jj >= O.startCol(ii) && jj < O.endCol(ii))
      {
        if (off_j == -1 && verticalBoundary.get(i, j - 1))
          continue;
        if (off_j == +1 && verticalBoundary.get(i, j))
          continue;
        if (off_i == -1 && (i == startRow || horizontalBoundary.get(i - 1, j)))
          continue;
        if (off_i == +1 && horizontalBoundary.get(i, j))
          continue;

        allNei[nNei] = T.get(ii, jj);
        sumNeighbor += allNei[nNei];
        nNei++;
      }
    }
  }
  allNei[nNei] = -1;
  return 0;
};

/*
This likelihood is used for p(ti|oi,Ni). Without unrelated T, Z.
It also discard o! directly. So it can only used to compare the relative likelihood of t 
at a position for a given fixed o. 
*/
double Estllike1(const double &o, const double &t, const double allNei[], const int &nNei, double bias)
{
  double variance = option_sigma;
  double lnSqXp1Xn1 = 0;
  for (int i = 0; i < nNei; i++)
  {
    lnSqXp1Xn1 += pow(log((t + 1) / (allNei[i] + 1)), 2);
  }
  return o * log(t * bias + 0.0001) - t * bias - lnSqXp1Xn1 / variance;
}

double greedyOptimizeT(const double &o, double t, const double allNei[], const int &nNei, double bias)
{
  double t_new;
  double dfdt;
  double lnXp1Xn1 = 0, lnSqXp1Xn1 = 0;
  double variance = option_sigma;
  double lr = 0.01;
  double f, newf;
  double maxIter = 10, iter = 0;
  for (int i = 0; i < nNei; i++)
  {
    lnXp1Xn1 += log((t + 1) / (allNei[i] + 1));
    lnSqXp1Xn1 += pow(log((t + 1) / (allNei[i] + 1)), 2);
  }
  dfdt = pow(bias * t, o - 1) * exp(-bias * t - lnSqXp1Xn1 / variance) * (o + pow(t, o) * (-bias - 2 * lnXp1Xn1 / (variance * (t + 1))));
  t_new = t + dfdt * lr;
  while (Estllike1(o, t_new, allNei, nNei, bias) < Estllike1(o, t, allNei, nNei, bias))
  {
    lr /= 5;
    t_new = t + dfdt * lr;
  }
  iter = 0;
  while (iter < maxIter && Estllike1(o, t_new, allNei, nNei, bias) > Estllike1(o, t, allNei, nNei, bias))
  {
    t = t_new;
    dfdt = pow(bias * t, o - 1) * exp(-bias * t - lnSqXp1Xn1 / variance) * (o + pow(t, o) * (-bias - 2 * lnXp1Xn1 / (variance * (t + 1))));
    t_new = t + dfdt * lr;
    iter++;
  }
  t = fmax(t, 0.0);
  return t;
}

double evaluateLikelihood(const float &currentT, const int &obs, const float &medNei, const float &bias)
{
  float transformedMedNeighbor;

  transformedMedNeighbor = log(medNei + 1);
  // no +1 here?
  float transformedCurrentTbias;
  transformedCurrentTbias = log((currentT + 0.00001) * bias);

  float localSigma = max(option_minSigma, sqrt(transformedMedNeighbor * option_sigmaMultiplier));
  float NORM_RES_over_sigma = NORM_RES / localSigma;

  float correctedp = currentT * bias;

  // log of Poisson distribution
  float loglike1 = -correctedp + obs * (transformedCurrentTbias)-logKFact[obs];
  return loglike1;

  float transformedp;
  transformedp = log(currentT + 1);

  float loglike2 = 0;
  /*int iii=(transformedp-transformedMedNeighbor)*NORM_RES_over_sigma + fullRange;
  if(iii<0) iii=0;
  else if(iii>floatFullRange) iii=floatFullRange;
  loglike2= precomputedNormLogPDF[iii]-transformedp;
*/
  float normalLogPDF = log((0.398942280401432678 / localSigma) * exp(-0.5 * pow((transformedp - transformedMedNeighbor) / localSigma, 2)));
  loglike2 = normalLogPDF - transformedp;
  return loglike1 + loglike2;
}

double evaluateFullLikelihood()
{
  double like = 0, l;
  double lModel = 0;
  double lObsGivenModel = 1;
  float t, medianNeighbor;
  int o;
  double biases;
  double variance = option_sigma;
  for (int i = O.startRow(); i < O.endRow(); i++)
  {
    for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
    {
      t = T.get(i, j);
      o = O.get(i, j);
      biases = bias[i] * bias[j];
      lObsGivenModel += -biases * t - logKFact[o] + o * log(biases * t + 0.00001);

      // lObsGivenModel += log(pow(biases * t, o) * exp(-biases * t) / KFact[o]);
      // cout<<pow(biases*t,o)<<" * "<<exp(-biases * t)<<" / "<<KFact[o]<<" "<<o<<endl;
      // // cout<<lObsGivenModel<<endl;
      // if(lObsGivenModel<-1e30)
      // exit(0);

      if (j + 1 < O.endCol(i) && !verticalBoundary.get(i, j))
        lModel -= pow(log(1 + T.get(i, j + 1)) - log(1 + t), 2) / variance;
      if (i + 1 < O.endRow())
      {
        if (j > O.startCol(i) + 1 && !horizontalBoundary.get(i, j))
          lModel -= pow(log(1 + T.get(i + 1, j)) - log(1 + t), 2) / variance;
        if (j + 1 < O.endCol(i + 1))
          if ((!horizontalBoundary.get(i, j) && !verticalBoundary.get(i + 1, j)) ||
              (!verticalBoundary.get(i, j) && !horizontalBoundary.get(i, j + 1)))
            lModel -= pow(log(1 + T.get(i + 1, j + 1)) - log(1 + t), 2) / variance;
        if (j - 1 > O.startCol(i + 1))
          if ((!verticalBoundary.get(i, j - 1) && !horizontalBoundary.get(i, j - 1)) ||
              (!horizontalBoundary.get(i, j) && !verticalBoundary.get(i + 1, j - 1)))
            lModel -= pow(log(1 + T.get(i + 1, j - 1)) - log(1 + t), 2) / variance;
      }
    }
  }
  return lObsGivenModel + lModel;
  // return like;
}
