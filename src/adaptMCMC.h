/*
 * adaptMCMC.h
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#ifndef adaptMCMC_H_
#define adaptMCMC_H_
#include "contactMap.h"
#include <queue>
#include <utility>
#include <map>

class adaptMCMC
{
public:
  int initialize(char *ofile);
  int burnIn();
  void outputSparseMatrixGZ(char *fn);
  int mainIter();

private:
  void outputCurrentSparseMatrixGZ(const char *fn);
  bool mainRuns;
  bool adaptive;
  UpperDiag<float> *sumOfTMat, *sumOfLogTMat, *sumOfSqTMat, *mContactMap11, *mContactMap12, *mContactMap13,
      *mContactMap21, *mContactMap22, *mContactMap23; // list of contactmaps to store data
  UpperDiag<float> *mCurrentTPtr, *mTPtr, *mT2Ptr,
      *mNext2Ptr, *mNextPtr, *mPrePtr, *mPre2Ptr, *mTmpPtr; //a list of pointer used to better working on data
  UpperDiag<char> *hasChanged;
  UpperDiag<unsigned char> *mAcceptances, *mAcceptances2, *mCurrentAcceptances; //,*chainB_Accept, *mCurrentChain_Accept;

  UpperDiag<float> *mSigmas, *mSigmas2, *mCurrentSigmas, *pairwisePotentialVariance; //scaleA, *scaleB, *mCurrentScale;

  UpperDiag<unsigned char> *perBaseChanges; // used in adaptMCMC in main runs
  unsigned long nChanges;
  float Adiff[9], Bdiff[9], Crossdiff[9];
  float getRegionalVar(int i,int j,int w);
  int diff;
  int It;
  int rep;
  int ppp;
  int chain;
  int nSamples, Changes;
  bool converaged, discover;
  int rowStep, colStep;
  char *ofile;
  double batchSize;
  std::queue<std::pair<int, int>> job;
  static void *blockMCMC(void *object)
  {
    reinterpret_cast<adaptMCMC *>(object)->blockMCMC();
    return 0;
  }

  static void *blockAdaptMCMC(void *object)
  {
    reinterpret_cast<adaptMCMC *>(object)->blockAdaptMCMC();
    return 0;
  }

  float mGetMedianNeighbor(const int &i, const int &j) const;

  int evaluateFullLikelihood(double& lPrior,double& lData);
  int getNeighbors(const int &i, const int &j, float allNei[]) const;

  double mhRatio(const int &ii, const int &jj, double &newVal);
  //   float saveAndTest();

  float mGetMedianNeighbor(const int &i, const int &j, const int &ip, const int &jp, const double &newVal, const bool &boundaries) const;

  int blockAdaptMCMC();

  int blockMCMC();
  int blockConverged();
  float sumChanges;
  int xxx;
  int OffSet;
  int rowOffsets[9] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
  int colOffsets[9] = {0, 1, 2, 0, 1, 2, 0, 1, 2};

  struct chainDiff
  {
    double diff[9][3]; // AA, BB, AB difference
    int counts;
  } mChainDiff;

  std::map<std::pair<int, int>, chainDiff> mBlockedChainDiff;
  double acceptLowBound, acceptHighBound;
  double *localPotentialSigma, *localPotential;

};

#endif /* SRC_MCMC_H_ */
