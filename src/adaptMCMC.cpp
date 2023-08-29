/*
 * MCMC.cpp
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */
#define MIN(a, b) ((a) > (b)) ? (b) : (a)
#define MAX(a, b) ((a) > (b)) ? (a) : (b)

#include "adaptMCMC.h"
#include "utils.h"
#include "readOptions.h"
#include "math.h"
#include "optimize.h"
#include "contactMap.h"
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <random>
#include "opt_median.h"
#include "poissonRegression.h"
#include "curvefit.h"
#include "bspline.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
pthread_mutex_t adaptmcmc_mutex;
thread_local bool th_boundaries = false;

float adaptMCMC::getRegionalVar(int i, int j, int w)
{
  // cout<<i<<" "<<j<<endl;
  double variance;
  double sumOfVal = 0, sumOfSqVal = 0;
  int n = 0;

  int ii, jj, off_j;
  int startRow = O.startRow(), endRow = O.endRow();
  double t, tNei;

  for (int off_i = -w; off_i <= w; off_i++)
  {
    for (off_j = -w; off_j <= w; off_j++)
    {
      ii = i + off_i;
      jj = j + off_j;
      // cout<<ii<<" "<<jj<<endl;
      // cout<<"startRow="<<startRow<<" endRow="<<endRow<<" O.startCol(ii)="<<O.startCol(ii)<<" O.endCol(ii)="<<O.endCol(ii)<<endl;
      if (ii < jj && ii >= startRow && ii < endRow && jj >= O.startCol(ii) && jj < O.endCol(ii))
      {
        // cout<<" -> "<<ii<<" "<<jj<<endl;

        t = log(T.get(ii, jj) + 1);
        if (ii - 1 >= startRow && jj - 1 >= O.startCol(ii - 1))
        {
          // cout<<"a\n"<<flush;
          // cout<<" -->"<<ii - 1<<","<< jj - 1<<endl;
          tNei = log(T.get(ii - 1, jj - 1) + 1);

          sumOfVal += fabs(tNei - t);
          sumOfSqVal += fabs((tNei - t) * (tNei - t));
          n++;
        }
        if (ii - 1 >= startRow && jj >= O.startCol(ii - 1) && jj < O.endCol(ii - 1))
        {
          // cout<<"b\n"<<flush;
          // cout<<" -->"<<ii - 1<<","<<jj<<endl;
          tNei = log(T.get(ii - 1, jj) + 1);
          sumOfVal += fabs(tNei - t);
          sumOfSqVal += fabs((tNei - t) * (tNei - t));
          n++;
        }
        if (jj - 1 >= O.startCol(ii) && ii < jj - 1)
        {
          // cout<<"c\n"<<flush;
          // cout<<" -->"<<ii<<","<<jj-1<<endl;
          tNei = log(T.get(ii, jj - 1) + 1);
          sumOfVal += fabs(tNei - t);
          sumOfSqVal += fabs((tNei - t) * (tNei - t));
          n++;
        }
        // cout<<exp(t)-1<<" tnei* "<<exp(tNei)-1<<endl;
      }
    }
  }
  // return 0.1;
  if (n == 0)
    return 1;
  // cout<<sumOfSqVal<<" "<<sumOfVal<<" "<<n<<endl;
  // else
  variance = (sumOfSqVal / n - (sumOfVal / n) * (sumOfVal / n)) * n / (n - 1);
  variance=fmax(0.1,variance);
  return variance;
}
int adaptMCMC::initialize(char *ofile)
{
  this->batchSize = 1;
  this->ofile = ofile;
  this->adaptive = true;
  double val;
  nSamples = 0;
  rowStep = stepSize;
  colStep = stepSize;
  // setBoundaries();
  // fprintf(stderr, "Done setting up boundaries\n");
  perBaseChanges = new UpperDiag<unsigned char>(option_firstRow,
                                                option_firstCol, option_lastRow,
                                                option_lastCol);

  //various contact map used for store intemidiate contacts in MCMC
  mContactMap11 = new UpperDiag<float>(option_firstRow, option_firstCol,
                                       option_lastRow, option_lastCol);
  mContactMap12 = new UpperDiag<float>(option_firstRow, option_firstCol,
                                       option_lastRow, option_lastCol);

  mContactMap21 = new UpperDiag<float>(option_firstRow, option_firstCol,
                                       option_lastRow, option_lastCol);
  mContactMap22 = new UpperDiag<float>(option_firstRow, option_firstCol,
                                       option_lastRow, option_lastCol);

  mSigmas = new UpperDiag<float>(option_firstRow, option_firstCol,
                                 option_lastRow, option_lastCol);
  mSigmas2 = new UpperDiag<float>(option_firstRow, option_firstCol,
                                  option_lastRow, option_lastCol);

  mAcceptances = new UpperDiag<unsigned char>(option_firstRow, option_firstCol,
                                              option_lastRow, option_lastCol);
  mAcceptances2 = new UpperDiag<unsigned char>(option_firstRow, option_firstCol,
                                               option_lastRow, option_lastCol);

  pairwisePotentialVariance = new UpperDiag<float>(option_firstRow, option_firstCol,
                                                   option_lastRow, option_lastCol);

  // initializing all Chains with rand.
  double sumOfT, sumOfC1, sumOfC2;
  sumOfT = 0;
  sumOfC1 = 0;
  sumOfC2 = 0;

  int w = option_sigma;
  double minPairwisePotentialVariance = 1e32;
  double currentPotentialVariance;
  for (int i = O.startRow(); i < O.endRow(); i++)
  {
    for (int j = O.startCol(i); j < O.endCol(i); j++)
    {
      currentPotentialVariance = getRegionalVar(i, j, w);
      pairwisePotentialVariance->set(i, j, currentPotentialVariance);
      if (currentPotentialVariance > 0 && currentPotentialVariance < minPairwisePotentialVariance)
        minPairwisePotentialVariance = currentPotentialVariance;
      // cout<<"pairwiseVariance: "<<i<<" "<<" "<<j<<" "<<pairwisePotentialVariance->get(i,j)<<endl;
    }
  }

  for (int i = O.startRow(); i < O.endRow(); i++)
  {
    for (int j = O.startCol(i); j < O.endCol(i); j++)
    {
      currentPotentialVariance = pairwisePotentialVariance->get(i, j);
      if (currentPotentialVariance < minPairwisePotentialVariance)
        pairwisePotentialVariance->set(i, j, minPairwisePotentialVariance);
      // cout << "pairwiseVariance: " << i << " "
      //      << " " << j << " " << pairwisePotentialVariance->get(i, j) << endl;
    }
  }
  // exit(0);

  for (int i = O.startRow(); i < O.endRow(); i++)
  {
    for (int j = O.startCol(i); j < O.endCol(i); j++)
    {
      mAcceptances->set(i, j, 0);
      mSigmas->set(i, j, max(0.01, T.get(i, j)));
      mSigmas->set(i, j, 1);
      mAcceptances2->set(i, j, 0);
      mSigmas2->set(i, j, max(0.01, T.get(i, j)));
      mSigmas2->set(i, j, 1);
      if (i != j)
      {
        sumOfT += T.get(i, j);
        val = drand48() * 1;
        mContactMap12->set(i, j, val);
        sumOfC1 += val;
        val = drand48() * 1;
        mContactMap22->set(i, j, val);
        sumOfC2 += val;
      }
      else
      {
        val = T.get(i, j);
        mContactMap12->set(i, j, val); //trick
        mContactMap22->set(i, j, val);
      }
    }
  }

  // //this is used to reset values near boundary under the bandSize execution;
  // //otherwise the algorithm may not converage due to it does not updating values outside the band
  // if (bandSize < O.endRow() - O.startRow())
  // {
  //   for (int i = O.startRow(); i < O.endRow(); i++)
  //   {
  //     for (int j = max(O.startCol(i), O.endCol(i) - colStep - 1); j < O.endCol(i); j++)
  //     {
  //       val = T.get(i, j);
  //       mContactMap12->set(i, j, val);
  //       mContactMap22->set(i, j, val);
  //     }
  //   }
  // }

  mPrePtr = mContactMap11;  //used to store previous contact map [kth before now] in burn-in
  mPre2Ptr = mContactMap21; //used to store previous contact map [kth before now] in burn-in

  mTPtr = mContactMap12;
  mT2Ptr = mContactMap22;

  sumOfTMat = new UpperDiag<float>(option_firstRow, option_firstCol, option_lastRow,
                                   option_lastCol);

  sumOfLogTMat = new UpperDiag<float>(option_firstRow, option_firstCol, option_lastRow,
                                      option_lastCol);

  sumOfSqTMat = new UpperDiag<float>(option_firstRow, option_firstCol, option_lastRow,
                                     option_lastCol);

  hasChanged = new UpperDiag<char>(option_firstRow, option_firstCol,
                                   option_lastRow, option_lastCol);

  nChanges = 1;
  diff = 0;
  converaged = false;

  mCurrentTPtr = mTPtr;

  //calculate local potential based on poisson regression
  localPotential = new UpperDiag<float>(option_firstRow, option_firstCol,
                                        option_lastRow, option_lastCol);

  double *breakpts = new double[5]; //cubic bspline with 7 ncoeffs (df=7-1)
  breakpts[0] = 0;
  breakpts[4] = breakpts[0];
  for (int i = O.startRow(); i < O.endRow(); i++)
    if (log(O.endCol(i) - i) > breakpts[4])
      breakpts[4] = log(O.endCol(i) - i);
  breakpts[1] = (breakpts[4] - breakpts[0]) * 0.25 + breakpts[0];
  breakpts[2] = (breakpts[4] - breakpts[0]) * 0.25 + breakpts[1];
  breakpts[3] = (breakpts[4] - breakpts[0]) * 0.25 + breakpts[2];

  cubicBspline cbs(7, breakpts);
  // exit(0);
  int numOfUsedPairs = 0;
  for (int i = O.startRow(); i < O.endRow(); i++)
    for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
      numOfUsedPairs++;
  double **X, *y;
  y = new double[max(510000, numOfUsedPairs)];
  X = new double *[max(510000, numOfUsedPairs)];
  for (int i = 0; i < max(510000, numOfUsedPairs); i++)
    X[i] = new double[8];
  cout << "numOfUsedPairs " << max(510000, numOfUsedPairs) << endl;

  int idx = 0;

  for (int i = O.startRow(); i < O.endRow(); i++)
  {
    for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
    {
      if (numOfUsedPairs > 500000 && drand48() > 500000.0 / numOfUsedPairs)
        continue;
      y[idx] = O.get(i, j);
      cbs.get_xi(log(j - i), &X[idx][0]);
      X[idx][7] = log(bias[i] * bias[j]);
      idx++;
    }
  }
  cout << "number of instance= " << idx << endl;
  poissonGLM pglm(8, idx, 1);
  cout << "start fitting\n";

  // pglm.fit(X, y);
  pglm.lbfgsfit(X, y);
  // exit(0);
  // cout << "predict \n";
  // for (int ttt = 0; ttt < 8; ttt++)
  //   cout << "X[0][i] " << X[100000][ttt] << endl;
  // cout << "===\n";
  // cout << pglm.predict(X[100000]);
  // exit(0);
  //refine null glm model
  double **X2, *y2;
  int idx2 = 0;
  y2 = new double[idx];
  X2 = new double *[idx];
  for (int i = 0; i < idx; i++)
    X2[i] = new double[8];
  double pval, mu;
  for (int i = 0; i < idx; i++)
  {
    mu = pglm.predict(X[i]);
    pval = gsl_cdf_poisson_Q(y[i], mu) + gsl_ran_poisson_pdf(y[i], mu);
    if (pval > 0.025)
    {
      y2[idx2] = y[i];
      for (int j = 0; j < 8; j++)
        X2[idx2][j] = X[i][j];
      idx2++;
    }
  }

  cout << "idx " << idx << endl;
  cout << "idx2 " << idx2 << endl;
  poissonGLM pglm2(8, idx2, 1);
  cout << "start refinement fitting\n";
  pglm2.lbfgsfit(X2, y2);

  delete[] y;
  delete[] y2;

  for (int i = 0; i < idx; i++)
  {
    delete[] X[i];
    delete[] X2[i];
  }
  delete[] X;
  delete[] X2;

  //. refine null glm nodel

  double x[8];
  for (int i = O.startRow(); i < O.endRow(); i++)
    for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
    {
      cout << "glm " << i << " " << j << " ";
      cbs.get_xi(log(j - i), x);
      x[7] = 0; //log(bias[i]*bias[j]);

      // x[0]=log(j-i);
      // x[1]=0;

      localPotential->set(i, j, pglm.predict(x));
      cout << pglm2.predict(x) << endl
           << flush;
    }

  //approx diagwise variance with normalized contact map(i.e. O/bias)
  double diagTsum, diagTTsum, logNormO;
  double *genomicdistance, *approxVar;
  genomicdistance = new double[max(O.endRow(), O.endCol(O.endRow() - 1))];
  approxVar = new double[max(O.endRow(), O.endCol(O.endRow() - 1))];
  localPotentialSigma = new double[max(O.endRow(), O.endCol(O.endRow() - 1))];
  int maxGenomicDistance = 0;
  int numOfValForVar = 0;
  int minVarNum = 20;
  idx = 0;
  double varlogNormO;
  for (int i = O.startRow(); i < O.endRow(); i++)
    for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
      maxGenomicDistance = max(j - i, maxGenomicDistance);

  for (int gd = 1; gd < maxGenomicDistance; gd++)
  {
    numOfValForVar = 0;
    diagTsum = 0;
    diagTTsum = 0;
    for (int i = O.startRow(); i < O.endRow(); i++)
      if (i + gd < O.endCol(i))
      {
        logNormO = log(O.get(i, i + gd) / (bias[i] * bias[i + gd]) + 1);
        // logNormO = (O.get(i, i + gd) / (bias[i] * bias[i + gd]));
        diagTTsum += logNormO * logNormO;
        diagTsum += logNormO;
        numOfValForVar++;
      }

    if (numOfValForVar > minVarNum)
    {
      varlogNormO = diagTTsum / numOfValForVar - pow(diagTsum / numOfValForVar, 2);
      varlogNormO *= numOfValForVar / (numOfValForVar - 1);
      if (varlogNormO > 0)
      {
        approxVar[idx] = varlogNormO;
        genomicdistance[idx] = gd;
        idx++;
        cout << "ttt " << gd << " " << varlogNormO << endl;
      }
    }
  }
  // return 0;
  double expparam_a, expparam_b, expparam_c;
  expparam_a = 0;
  expparam_b = 0;
  expparam_c = 0;

  // fitExp_bfgs(idx, genomicdistance, approxVar, expparam_a, expparam_b, expparam_c);
  fitExp(idx, genomicdistance, approxVar, expparam_a, expparam_b, expparam_c);
  cerr << "y=" << expparam_a << "*np.exp(" << expparam_b << "*x)+" << expparam_c << endl;
  // exit(0);
  delete[] genomicdistance;
  delete[] approxVar;
  // cout << expparam_a << " " << expparam_b << " " << expparam_c << endl;

  for (int i = 0; i < max(O.endRow(), O.endCol(O.endRow() - 1)); i++)
  {
    // cout << "localSigmaSq " << i << " ";
    localPotentialSigma[i] = predExp(i, expparam_a, expparam_b, expparam_c);
    localPotentialSigma[i] = fmax(1e-5, localPotentialSigma[i]);
    // cout << localPotentialSigma[i] << endl
    //  << flush;
  }
  cout << "finished localPotentialSigma\n";

  return 0;
}

int adaptMCMC::burnIn()
{
  cout << "start burn-in\n";
  int obs;
  double val, newVal;
  double ll1, ll2;

  mainRuns = false;
  double lPrior, lData;
  this->mChainDiff.counts = 0;
  memset(this->mChainDiff.diff, 0, sizeof(this->mChainDiff.diff));
  int startRow = O.startRow();
  int endRow = O.endRow();

  if (bandSize < O.endRow() - O.startRow())
  {
    for (int i = startRow; i < endRow; i += rowStep)
    {
      for (int j = O.startCol(i); j < O.endCol(min(i + colStep - 1, endRow)); j += colStep)
      {
        this->mBlockedChainDiff[make_pair(i, j)] = this->mChainDiff;
      }
    }
  }
  else
  {
    for (int i = startRow; i < endRow; i += rowStep)
    {
      for (int j = O.startCol(i); j < O.endCol(i); j += colStep)
      {
        this->mBlockedChainDiff[make_pair(i, j)] = this->mChainDiff;
      }
    }
  }
  int tids[n_threads];
  pthread_t threads[n_threads];
  pthread_mutex_init(&adaptmcmc_mutex, NULL);
  for (rep = 1; true; rep++)
  {
    // cout << "burn-in iter: " << rep << endl;
    if (rep == 1300)
      return 0;
    // this->adaptive=false;
    for (this->chain = 1; this->chain <= 2; this->chain++)
    {
      Changes = 0;
      sumChanges = 0;
      if (chain == 1)
      {
        mCurrentTPtr = mTPtr;
        mCurrentAcceptances = mAcceptances;
        mCurrentSigmas = mSigmas;
      }
      else
      {
        mCurrentTPtr = mT2Ptr;
        mCurrentAcceptances = mAcceptances2;
        mCurrentSigmas = mSigmas2;
      }

      startRow = O.startRow();
      endRow = O.endRow();
      for (int i = startRow; i < endRow; i += rowStep)
      {

        if (bandSize > O.endRow() - O.startRow())
        {
          for (int j = O.startCol(i); j < O.endCol(i); j += colStep)
          {
            job.push(make_pair(i, j));
          }
        }
        else
        {
          for (int j = O.startCol(i); j < O.endCol(min(i + colStep - 1, endRow)); j += colStep)
          {
            job.push(make_pair(i, j));
          }
        }
      }

      for (int i = 0; i < n_threads; i++)
      {
        tids[i] = i;
        pthread_create(&threads[i], NULL, blockMCMC, this);
      }

      for (int i = 0; i < n_threads; i++)
      {
        pthread_join(threads[i], NULL);
      }

      if (chain == 1)
      {
        if (rep % 50 == 1)
        {
          this->evaluateFullLikelihood(lPrior, lData);
          cout << "burnIn it" << rep << ",chain1, Likelihood: " << lPrior + lData << " (" << lPrior << " , " << lData << ") changes=" << Changes << "\tsumChanges: " << sumChanges;
        }
      }
      else
      {
        if (rep % 50 == 1)
        {
          this->evaluateFullLikelihood(lPrior, lData);
          cout << " ;chain2 changes=" << Changes << " , sumChanges " << sumChanges
               << " Likelihood: " << lPrior + lData << " (" << lPrior << " , " << lData << ")" << endl;
        }
      }
    } //end of looping chain 1, 2

    //this is used for adaptMCMC
    if (this->adaptive && rep % 50 == 0 && rep > 50)
    {
      for (int i = O.startRow(); i < O.endRow(); i++)
      {
        for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
        {

          if (mAcceptances->get(i, j) - 0 < 14 && mSigmas->get(i, j) > 0.01)
          {

            mSigmas->set(i, j, this->mSigmas->get(i, j) * 0.5);
          }

          else if (mAcceptances->get(i, j) - 0 > 30 && mSigmas->get(i, j) < 10000)
          {

            mSigmas->set(i, j, this->mSigmas->get(i, j) * 1.5);
          }
          mAcceptances->set(i, j, 0);

          if (mAcceptances2->get(i, j) - 0 < 14 && mSigmas2->get(i, j) > 0.01)
          {

            mSigmas2->set(i, j, mSigmas2->get(i, j) * 0.5);
          }
          else if (mAcceptances2->get(i, j) - 0 > 30 && mSigmas2->get(i, j) < 10000)
          {

            mSigmas2->set(i, j, this->mSigmas2->get(i, j) * 1.5);
          }
          mAcceptances2->set(i, j, 0);
        }
      }
    } //end of update sigma in adaptMCMC

    //convergence test after each k iterations
    if (rep == k)
    {
      cout << "Burn-in iter: " << rep << endl;
      for (int i = O.startRow(); i < O.endRow(); i++)
      {
        for (int j = O.startCol(i); j < O.endCol(i); j++)
        {
          mPrePtr->set(i, j, mTPtr->get(i, j));
          mPre2Ptr->set(i, j, mT2Ptr->get(i, j));
        }
      }
      diff = 0;
    }
    else if (rep > 0 && rep % k == 0)
    {
      cout << "Burn-in iter: " << rep << endl;
      if (blockConverged())
      {
        cout << "burn-in done\n";
        mCurrentTPtr = mTPtr;
        return 0;
      }
      for (int i = O.startRow(); i < O.endRow(); i++)
      {
        for (int j = O.startCol(i); j < O.endCol(i); j++)
        {
          mPrePtr->set(i, j, mTPtr->get(i, j));
          mPre2Ptr->set(i, j, mT2Ptr->get(i, j));
        }
      }
    }
  } // end of MCMC iteration loop
  return 0;
}

int adaptMCMC::mainIter()
{
  mainRuns = true;
  this->chain = 1;
  mCurrentTPtr = mTPtr;
  mCurrentAcceptances = mAcceptances;
  mCurrentSigmas = mSigmas;
  double lPrior, lData;
  int startRow = O.startRow(), endRow = O.endRow();
  this->evaluateFullLikelihood(lPrior, lData);
  cerr << "MCMC iter: " << 0 << " Likelihood: " << lPrior + lData << " (" << lPrior << " , " << lData << ") Changes:"
       << Changes;
  if (option_testMatrix[0])
    cerr << ", test_SSE: " << evaluateAgainstTest(mCurrentTPtr);
  cerr << endl;
  int tids[n_threads];
  double newVal;
  pthread_t threads[n_threads];
  pthread_mutex_init(&adaptmcmc_mutex, NULL);
  for (rep = 0; rep < option_mcmcMaxIter; rep++)
  {
    sumChanges = 0;
    this->It = rep;
    Changes = 0;

    // for (int i = O.startRow(); i < O.endRow(); i += rowStep)
    // {
    //   for (int j = O.startCol(i); j < O.endCol(i); j += colStep)
    //   {
    //     job.push(make_pair(i, j));
    //   }
    // }

    for (int i = startRow; i < endRow; i += rowStep)
    {

      if (bandSize > O.endRow() - O.startRow())
      {
        for (int j = O.startCol(i); j < O.endCol(i); j += colStep)
        {
          job.push(make_pair(i, j));
        }
      }
      else
      {
        for (int j = O.startCol(i); j < O.endCol(min(i + colStep - 1, endRow)); j += colStep)
        {
          job.push(make_pair(i, j));
        }
      }
    }

    for (int i = 0; i < n_threads; i++)
    {
      tids[i] = i;
      pthread_create(&threads[i], NULL, blockMCMC, this);
    }

    for (int i = 0; i < n_threads; i++)
    {
      pthread_join(threads[i], NULL);
    }

    if (rep % 10 == 0)
    {
      nSamples++;
      for (int i = O.startRow(); i < O.endRow(); i++)
      {
        for (int j = O.startCol(i); j < O.endCol(i); j++)
        {
          newVal = mCurrentTPtr->get(i, j);
          this->sumOfLogTMat->add(i, j, log(newVal));
          this->sumOfTMat->add(i, j, newVal);
          this->sumOfSqTMat->add(i, j, newVal * newVal);
        }
      }
      if (outputSamples)
      {
        //here we print collected samples
        std::stringstream ss;
        ss << this->ofile << ".rep_" << rep << ".sparse.gz";
        this->outputCurrentSparseMatrixGZ(ss.str().c_str());
      }
    }
    /*  adaptive is not allowed after the chain converged
    if (rep % 50 == 0 && rep >= 50)
    {
      for (int i = O.startRow(); i < O.endRow(); i++)
      {

        for (int j = O.startCol(i); j < O.endCol(i); j++)
        {

          if (mCurrentAcceptances->get(i, j) - 0 < 14 && mCurrentSigmas->get(i, j) > 0.01)
          {

            mCurrentSigmas->set(i, j, this->mCurrentSigmas->get(i, j) * 0.5);
          }
          else if (mCurrentAcceptances->get(i, j) - 0 > 30 && mCurrentSigmas->get(i, j) < 10000)
          {

            mCurrentSigmas->set(i, j, this->mCurrentSigmas->get(i, j) * 1.5);
          }
          mCurrentAcceptances->set(i, j, 0);
        }
      }
    }
    */
    if (rep % 100 == 0)
    {
      this->evaluateFullLikelihood(lPrior, lData);
      cerr << "MCMC iter: " << rep + 1 << " Likelihood: " << lPrior + lData << " (" << lPrior << " , " << lData << ") Changes:"
           << Changes << ", sumChanges: " << sumChanges;
      if (option_testMatrix[0])
        cerr << ", test_SSE: " << evaluateAgainstTest(mCurrentTPtr);
      cerr << endl;
    }
  } // end of MCMC iteration loop

  for (int i = O.startRow(); i < O.endRow(); i++)
  {
    for (int j = O.startCol(i); j < O.endCol(i); j++)
    {
      T.set(i, j, mCurrentTPtr->get(i, j));
    }
  } // copy final result to T matrix
  return 0;
}

// MCMC in each block
int adaptMCMC::blockMCMC()
{
  double val, newVal, ll1, ll2, medij;
  int obs;
  double randratio = 0;
  int startRow, startCol, endRow, endCol;
  bool changed;
  int Jstart, Jend;
  double rr;
  while (true)
  {

    if (job.empty())
      return 0;
    pthread_mutex_lock(&adaptmcmc_mutex);
    if (job.empty())
    {
      pthread_mutex_unlock(&adaptmcmc_mutex);
      return 0;
    }

    startRow = job.front().first;
    startCol = job.front().second;
    job.pop();
    pthread_mutex_unlock(&adaptmcmc_mutex);

    endRow = MIN(O.endRow(), startRow + this->rowStep + 1);
    std::vector<std::pair<int, int>> entries2DIdx;
    for (int i = startRow; i < endRow; i++)
    {
      Jstart = MAX(startCol, O.startCol(i));
      Jend = MIN(startCol + this->colStep + 1, O.endCol(i)); // jag-jag for band version
      for (int j = Jstart; j < Jend; j++)
      {
        entries2DIdx.push_back(make_pair(i, j));
      }
    }
    shuffle(entries2DIdx.begin(), entries2DIdx.end(), default_random_engine(0));

    int i, j;
    while (!entries2DIdx.empty())
    // for (int j = Jstart; j < Jend; j++)
    {
      i = entries2DIdx.back().first;
      j = entries2DIdx.back().second;
      entries2DIdx.pop_back();
      val = mCurrentTPtr->get(i, j);

      obs = O.get(i, j);
      if (i == j)
      {
        continue;
      }
      newVal = exp(log(val) + norm_Randn(0, mCurrentSigmas->get(i, j)));
      // if(i==200 && j==234)
      // cout<<mCurrentSigmas->get(i,j)<<endl;
      // newVal = exp(log(val) + norm_Randn(0, 1));

      rr = drand48();

      if (mhRatio(i, j, newVal) > rr)
      {
        Changes++;
        mCurrentTPtr->set(i, j, newVal);
        mCurrentAcceptances->increase(i, j);
        sumChanges += (newVal - val) * (newVal - val);
      }
    }
  }
  return 0;
}

void adaptMCMC::outputSparseMatrixGZ(char *fn)
{
  double mean, variance, k, theta, s;
  gzFile gzO = gzopen(fn, "w");
  gzprintf(gzO, "#nSamples= %d\n", nSamples);
  gzprintf(gzO, "chr1\tpos1\tchr2\tpos2\tobs\tbias\tmean\tvariance\ttheta\tk\n", nSamples);
  double eps = 1e-12;

  for (int i = T.startRow(); i < T.endRow(); i++)
  {

    for (int j = T.startCol(i); j < T.endCol(i); j++)
    {
      mean = sumOfTMat->get(i, j) / nSamples;
      variance = sumOfSqTMat->get(i, j) / nSamples - pow(mean, 2);

      s = log(mean) - sumOfLogTMat->get(i, j) / nSamples;
      k = (3 - s + sqrt(pow(s - 3, 2) + 24 * s)) / (12 * s + eps);
      theta = sumOfTMat->get(i, j) / ((k + eps) * nSamples);

      // if (!option_outputNormalized)
      // {
      //   mean = S->get(i, j) / nSamples * (bias[i] * bias[j]);
      //   variance = (SS->get(i, j) / (nSamples - 1) * (bias[i] * bias[j]) * (bias[i] * bias[j]) - pow(mean, 2));
      // }
      // else
      // {

      // }
      if (mean > 0)
        gzprintf(gzO, "%d\t%d\t%d\t%d\t%5.5lf\t%5.5lf\t%5.5lf\t%5.5lf\t%5.5lf\t%5.5lf\n", chr1, i, chr2, j, O.get(i, j), bias[i] * bias[j], mean, variance, theta, k);
    }
  }
  gzclose(gzO);
}

float adaptMCMC::mGetMedianNeighbor(const int &i, const int &j) const
{
  int nNei = 0;
  float sumNeighbor = 0;
  unsigned short int seed;
  float allNei[8];
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

        allNei[nNei] = exp(mCurrentTPtr->get(ii, jj));

        sumNeighbor += allNei[nNei];
        nNei++;
      }
    }
  }
  if (option_median)
  {

    switch (nNei)
    {
    case 0:
      return 0;
    case 1:
      return allNei[0];
    case 2:
      return (allNei[0] + allNei[1]) / 2;
    case 3:
      return opt_med3(allNei);
    case 4:
      return opt_med4(allNei);
    case 5:
      return opt_med5(allNei);
    case 6:
      return opt_med6(allNei);
    case 7:
      return opt_med7(allNei);
    case 8:
      return opt_med8(allNei);
    case 9:
      return opt_med9(allNei);
    }
  }
  return sumNeighbor / nNei;
}

int adaptMCMC::evaluateFullLikelihood(double &lPrior, double &lData)
{
  double biases;
  int nNei = 0;
  float sumNeighbor = 0;
  float *allNei;
  allNei = new float[9];

  double o;
  double t;
  double variance = option_sigma;
  lPrior = 0;
  lData = 0;

  for (int i = O.startRow(); i < O.endRow(); i++)
  {
    for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
    {
      t = mCurrentTPtr->get(i, j);
      o = O.get(i, j);
      biases = bias[i] * bias[j];
      getNeighbors(i, j, allNei);
      for (int ii = 0; ii < 8 && allNei[ii] > -1; ii++)
      {

        lPrior -= pow(log(allNei[ii] + 1) - log(t + 1), 2) / pairwisePotentialVariance->get(i, j); // variance; //fmax(0.1,variance*fmax(log(allNei[ii] + 1),log(oldVal + 1)));
      }

      lPrior -= pow(log(localPotential->get(i, j) + 1) - log(t + 1), 2) / localPotentialSigma[j - i];
      //poisson
      lData += o * log(biases * t * this->batchSize) - biases * t * this->batchSize;
      //end of poisson
    }
  }

  delete[] allNei;
  return 0;
}

int adaptMCMC::getNeighbors(const int &i, const int &j, float allNei[]) const
{
  int nNei = 0;
  float sumNeighbor = 0;
  int ii, jj, off_j;
  int startRow = O.startRow(), endRow = O.endRow();
  if (i + 1 == j)
  {
    if (i - 1 >= startRow)
      allNei[nNei++] = mCurrentTPtr->get(i - 1, j - 1);
    if (i + 1 < endRow)
      allNei[nNei++] = mCurrentTPtr->get(i + 1, j + 1);
    allNei[nNei] = -1;
    return 0;
  }
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
        // if (off_j == -1 && verticalBoundary.get(i, j - 1))
        //   continue;
        // if (off_j == +1 && verticalBoundary.get(i, j))
        //   continue;
        // if (off_i == -1 && (i == startRow || horizontalBoundary.get(i - 1, j)))
        //   continue;
        // if (off_i == +1 && horizontalBoundary.get(i, j))
        //   continue;

        // if(ii==i || jj==j){ //used for first order
        allNei[nNei] = mCurrentTPtr->get(ii, jj);
        sumNeighbor += allNei[nNei];
        nNei++;
        // }
      }
    }
  }
  allNei[nNei] = -1;
  return 0;
};

/***************
 * 
 * 
 * 
 * *************** */
//return p(new)/p(old);
double adaptMCMC::mhRatio(const int &i, const int &j, double &newVal)
{
  double llike1, llike2;
  double biases = bias[i] * bias[j];
  int nNei = 0;
  float sumNeighbor = 0;
  float *allNei;
  allNei = new float[9];
  getNeighbors(i, j, allNei);
  double oldVal = mCurrentTPtr->get(i, j);
  double o = O.get(i, j);
  double variance; // = option_sigma;
  llike1 = 0;
  llike2 = 0;
  for (int ii = 0; ii < 8 && allNei[ii] > -1; ii++)
  {
    // cout<<"variance "<<variance<<endl;
    llike1 -= pow(log(allNei[ii] + 1) - log(oldVal + 1), 2) / pairwisePotentialVariance->get(i, j); //variance; //fmax(0.1,variance*fmax(log(allNei[ii] + 1),log(oldVal + 1)));
    llike2 -= pow(log(allNei[ii] + 1) - log(newVal + 1), 2) / pairwisePotentialVariance->get(i, j); //variance; //fmax(0.1,variance*fmax(log(allNei[ii] + 1),log(newVal + 1)));
  }

  //poisson
  llike1 += o * log(biases * oldVal * this->batchSize) - biases * oldVal * this->batchSize;
  llike2 += o * log(biases * newVal * this->batchSize) - biases * newVal * this->batchSize;
  //end of poisson

  // // // //local potential
  // // // // gamma based local potential.. removed
  // // // llike1 += (piror_k[abs(j - i)] - 1) * log(oldVal) - oldVal / piror_theta[abs(j - i)];
  // // // llike2 += (piror_k[abs(j - i)] - 1) * log(newVal) - newVal / piror_theta[abs(j - i)];
  // // // //==================================

  llike1 -= pow(log(localPotential->get(i, j) + 1) - log(oldVal + 1), 2) / localPotentialSigma[j - i];
  llike2 -= pow(log(localPotential->get(i, j) + 1) - log(newVal + 1), 2) / localPotentialSigma[j - i];

  // llike1 -= pow(localPotential->get(i, j) - oldVal * this->batchSize, 2) / localPotentialSigma[j - i];
  // llike2 -= pow(localPotential->get(i, j) - newVal * this->batchSize, 2) / localPotentialSigma[j - i];

  delete[] allNei;
  return exp(llike2 - llike1) * newVal / oldVal;
}

int adaptMCMC::blockConverged()
{
  if (this->mBlockedChainDiff.empty())
    return 1;

  int startRow = O.startRow();
  int endRow = O.endRow();
  int startCol, Jend, Jstart;
  int startRow2, startCol2, endRow2, endCol2, Jend2, Jstart2;
  int OendCol;
  std::pair<int, int> key;
  chainDiff *cdiff;
  float Adiff_ = 0, Bdiff_ = 0, crossdiff_ = 0;
  float sumA, sumB, sumCross;
  float val11, val12, val21, val22;
  for (int i = startRow; i < endRow; i += this->rowStep)
  {
    if (bandSize < O.endRow() - O.startRow())
      OendCol = O.endCol(min(i + colStep - 1, endRow));
    else
      OendCol = O.endCol(i);
    for (int j = O.startCol(i); j < OendCol; j += this->colStep)
    {
      key = make_pair(i, j);
      if (this->mBlockedChainDiff.find(key) == mBlockedChainDiff.end())
      {
        continue;
      }
      cdiff = &(this->mBlockedChainDiff.find(key)->second);
      cout << i << "," << j << " ~ " << i + this->rowStep << ","
           << j + this->colStep << endl;
      //shift 1 position left of difference for incoorp new one
      if (cdiff->counts == 9)
      {
        for (int diffIdx = 1; diffIdx < cdiff->counts; diffIdx++)
        {
          cdiff->diff[diffIdx - 1][0] = cdiff->diff[diffIdx][0];
          cdiff->diff[diffIdx - 1][1] = cdiff->diff[diffIdx][1];
          cdiff->diff[diffIdx - 1][2] = cdiff->diff[diffIdx][2];
        }
        cdiff->counts--;
      }
      sumA = 0;
      sumB = 0;
      sumCross = 0;

      Adiff_ = 0;
      Bdiff_ = 0;
      crossdiff_ = 0;
      startRow2 = i;
      startCol2 = j;
      int totalEntries = 0;
      endRow2 = MIN(O.endRow(), startRow2 + this->rowStep);

      for (int ii = startRow2; ii < endRow2; ii++)
      {
        Jstart2 = MAX(startCol2, O.startCol(ii)) + 1;
        for (int jj = Jstart2;
             jj < min(startCol2 + this->colStep, O.endCol(ii)); jj++)
        {
          totalEntries++;
          val11 = mPrePtr->get(ii, jj);
          val12 = mTPtr->get(ii, jj);
          val21 = mPre2Ptr->get(ii, jj);
          val22 = mT2Ptr->get(ii, jj);

          Adiff_ += pow(val11 - val12, 2);
          Bdiff_ += pow(val21 - val22, 2);
          crossdiff_ += pow(val22 - val12, 2);
        }
      }

      cdiff->diff[cdiff->counts][0] = (Adiff_ / totalEntries);
      cdiff->diff[cdiff->counts][1] = (Bdiff_ / totalEntries);
      cdiff->diff[cdiff->counts][2] = (crossdiff_ / totalEntries);
      cdiff->counts += 1;
      // if (cdiff->counts >0)
      if (cdiff->counts == 9)
      {
        for (int diffIdx = 0; diffIdx < cdiff->counts; diffIdx++)
        {
          sumA += cdiff->diff[diffIdx][0];
          sumB += cdiff->diff[diffIdx][1];
          sumCross += cdiff->diff[diffIdx][2];
        }
        sumA /= cdiff->counts;
        sumB /= cdiff->counts;
        sumCross /= cdiff->counts;

        if (sumA == 0 && sumB == 0 && sumCross == 0)
        {
          cout << "iter: " << this->It << " " << i << "," << j << " ~ " << i + this->rowStep << ","
               << j + this->colStep;
          cout << " converaged: " << 0 << endl;
          this->mBlockedChainDiff.erase(key);
        }
        else if (fabs(sumCross - sumA) <= 0.1 * sumA)
        // else if (fabs(sumCross - sumA) <= -0.1 * sumA)
        {
          cout << "iter: " << this->It << " " << i << "," << j << " ~ " << i + this->rowStep << ","
               << j + this->colStep;
          cout << " converaged: "
               << fabs(sumCross - sumA) / sumA << endl;

          this->mBlockedChainDiff.erase(key);
        }
        else
        {
          cout << "iter: " << this->It << ":\t";
          cout << "inter- " << sumCross << "\t";
          cout << "intra- " << sumA << " , " << sumB << "\t";
          if (sumCross < 1e-3 && sumA < 1e-3 && sumB < 1e-3)
          {
            cout << " converaged\n";
            this->mBlockedChainDiff.erase(key);
          }
          else
            cout << " not converaged: " << fabs(sumCross - sumA) / sumA << endl;
        }
      }
    }
  }
  cout << "non-converaged blocks=" << this->mBlockedChainDiff.size() << endl;
  if (this->mBlockedChainDiff.empty())
    return 1;
  return 0;
}

float adaptMCMC::mGetMedianNeighbor(const int &i, const int &j, const int &ip, const int &jp,
                                    const double &newVal, const bool &boundaries) const
{
  int nNei = 0;
  float sumNeighbor = 0;
  float allNei[8];
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
      if (boundaries && ii < jj && ii >= startRow && ii < endRow && jj >= O.startCol(ii) && jj < O.endCol(ii))
      {
        if (off_j == -1 && verticalBoundary.get(i, j - 1))
          continue;
        if (off_j == +1 && verticalBoundary.get(i, j))
          continue;
        if (off_i == -1 && (i == startRow || horizontalBoundary.get(i - 1, j)))
          continue;
        if (off_i == +1 && horizontalBoundary.get(i, j))
          continue;
        if (ii == ip && jj == jp)
        {
          allNei[nNei] = newVal;
        }
        else
          allNei[nNei] = mCurrentTPtr->get(ii, jj);
        sumNeighbor += allNei[nNei];
        nNei++;
      }
      else if (!boundaries)
      {
        allNei[nNei] = mCurrentTPtr->get(ii, jj);
        sumNeighbor += allNei[nNei];
        nNei++;
      }
    }
  }

  if (option_median)
  {
    switch (nNei)
    {
    case 0:
      return 0;
    case 1:
      return allNei[0];
    case 2:
      return (allNei[0] + allNei[1]) / 2;
    case 3:
      return opt_med3(allNei);
    case 4:
      return opt_med4(allNei);
    case 5:
      return opt_med5(allNei);
    case 6:
      return opt_med6(allNei);
    case 7:
      return opt_med7(allNei);
    case 8:
      return opt_med8(allNei);
    case 9:
      return opt_med9(allNei);
    }
  }
  else
    return sumNeighbor / nNei;
}

void adaptMCMC::outputCurrentSparseMatrixGZ(const char *fn)
{
  float val;
  gzFile gzO = gzopen(fn, "w");
  // if (option_lastRow != -1 && option_lastCol != -1)
  //   gzprintf(gzO, "# %d %d %d %d\n", option_firstRow, option_lastRow,
  //            option_firstCol, option_lastCol);
  // else
  //   gzprintf(gzO, "# %d %d %d %d\n", option_firstRow, option_lastRow,
  //            option_firstCol, option_lastCol);

  for (int i = T.startRow(); i < T.endRow(); i++)
  {

    for (int j = T.startCol(i); j < T.endCol(i); j++)
    {
      if (i == j)
        continue;
      val = this->mCurrentTPtr->get(i, j);
      //if (val > option_minOutput)
      gzprintf(gzO, "%d %d %d %d %5.4lf\n",chr1, i,chr2, j,val);
    }
  }
  gzclose(gzO);
}
