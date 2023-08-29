/*
 * KDE.cpp
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#include "KDE.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "contactMap.h"
#include "readOptions.h"
#include "utils.h"
using namespace std;

pthread_mutex_t kde_mutex;
int nbStdDev = 3;

KDE::KDE()
{
	// TODO Auto-generated constructor stub
}

KDE::~KDE()
{
	// TODO Auto-generated destructor stub
}

int KDE::initialize()
{
	cumO = new UpperDiag<unsigned long>(option_firstRow, option_firstCol,
										option_lastRow, option_lastCol);
	nextNonEmpty = new UpperDiag<int>(option_firstRow, option_firstCol,
									  option_lastRow, option_lastCol);
	// build sparse matrix representation
	fprintf(stderr, "Setting up nextNon-empty\n");
	fflush(stderr);
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		int last = O.startCol(i);
		for (int j = O.startCol(i); j < O.endCol(i); j++)
		{
			nextNonEmpty->set(i, j, -1);
			if (O.get(i, j) != 0)
			{
				//			if (O.get(i, j) != 0 && i!=j) { // i!=j means skip self interactions
				for (int jj = last; jj < j; jj++)
				{
					nextNonEmpty->set(i, jj, j);
				}
				last = j;
			}
		}
	}

	//   initialize T matrix
	fprintf(stderr, "Setting up T\n");
	int sumObs = 0;

	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i); j < O.endCol(i); j++)
		{
			//		  if(i==j)// self interaction
			//		    continue;
			T.set(i, j, -1);
			sumObs += O.get(i, j);
			//      if (t) fprintf(stderr,"O[%d][%d]=%d\n",i,j,t);
		}
	}
	fprintf(stderr, "sumObs = %d\n", sumObs);

	mini_h = 0.3;
	maxi_h = option_kdeMaxBandwidth;

	if (option_kdeBandwidth != -1)
	{
		mini_h = option_kdeBandwidth;
		maxi_h = mini_h + 0.0001;
	}

	// precomputing OoverBias matrix
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i); j < O.endCol(i); j++)
		{
			//		  if (i==j)
			//		    continue;
			int x = O.get(i, j);
			if (x > 0)
				OoverBias.set(i, j, x / (bias[i] * bias[j]));
		}
	}

	multFactor = 1.3;
	lastIter = false;

	rowStep = 300;
	colStep = 300;

	return 0;
}

int KDE::runKDE()
{
	int tids[n_threads];
	pthread_t threads[n_threads];
	pthread_mutex_init(&kde_mutex, NULL);
	for (float h = mini_h; h < maxi_h; h *= multFactor)
	{
		fprintf(stderr, "h=%f\n", h);
		if (h * multFactor >= maxi_h)
			lastIter = true;
		nbStdDevTh = (int)(nbStdDev * h);
		// compute cumulative: cumO->set(i, j, above + left - aboveleft + O.get(i, j));
		fprintf(stderr, "Computing cumulative\n");
		if (h == mini_h)
		{
			//				computeCumulative();
			for (int i = O.startRow(); i < O.endRow(); i++)
			{
				for (int j = O.startCol(i); j < O.endCol(i); j++)
				{
					//				  if (i==j)
					//				    continue;
					unsigned long above =
						(i > O.startRow()) ? cumO->get(i - 1, j) : 0;

					unsigned long left = 0;
					if (j > O.startCol(i))
						left = cumO->get(i, j - 1);
					else
					{
						if (i > O.startRow() && j - 1 >= O.startCol(i - 1))
							left = cumO->get(i - 1, j - 1); //this is above left
					}

					unsigned long aboveleft = 0;
					if (i - 1 >= O.startRow() && j - 1 >= O.startCol(i - 1))
						aboveleft = cumO->get(i - 1, j - 1);

					cumO->set(i, j, above + left - aboveleft + O.get(i, j));
				}
			}
		}

		// 2D gaussian distribution
		fprintf(stderr, "Computing gaussian\n");
		gaussian = (float **)malloc(
			sizeof(float *) * (int)(nbStdDev * 2 * (maxi_h + 1) + 1));
		for (int i = 0; i < (int)(nbStdDev * 2 * (maxi_h + 1) + 1); i++)
			gaussian[i] = (float *)malloc(
				sizeof(float) * (int)((nbStdDev * 2 * (maxi_h + 1) + 1)));

		float s = 0;
		for (int i = -nbStdDevTh; i <= nbStdDevTh; i++)
		{
			for (int j = -nbStdDevTh; j <= nbStdDevTh; j++)
			{
				gaussian[i + nbStdDevTh][j + nbStdDevTh] = 0;
				normpdf(sqrt(i * i + j * j), 0, h);
				gaussian[i + nbStdDevTh][j + nbStdDevTh] = normpdf(
					sqrt(i * i + j * j), 0, h);
				s += gaussian[i + nbStdDevTh][j + nbStdDevTh];
			}
		}

		for (int i = -nbStdDevTh; i <= nbStdDevTh; i++)
		{
			for (int j = -nbStdDevTh; j <= nbStdDevTh; j++)
			{
				gaussian[i + nbStdDevTh][j + nbStdDevTh] /= s;
			}
		}

		nSet = 0;
		maxRadius2 = nbStdDevTh * nbStdDevTh;
		maxRadius = nbStdDevTh;

		minb = (int *)malloc(sizeof(int) * (maxRadius * 2 + 3));
		for (int a = -maxRadius; a <= maxRadius; a++)
		{
			minb[(int)(a + maxRadius)] = (int)(sqrt(maxRadius2 - a * a));
		}

		gaussianSum = 0;
		for (int a = -maxRadius; a <= maxRadius; a++)
		{
			int range = minb[(int)(a + maxRadius)];
			for (int b = -range; b <= range; b++)
			{
				gaussianSum += gaussian[a + nbStdDevTh][b + nbStdDevTh];
			}
		}

		// iterate over matrix to apply gaussian filter
		for (int i = O.startRow(); i < O.endRow(); i += rowStep)
		{
			for (int j = O.startCol(i); j < O.endCol(i); j += colStep)
			{
				job.push(make_pair(i, j));
			}
		}
		cout << "size(job):" << job.size() << endl;
		for (int i = 0; i < n_threads; i++)
		{
			tids[i] = i;
			pthread_create(&threads[i], NULL, applyGaussianFilter, this);
		}

		for (int i = 0; i < n_threads; i++)
		{
			pthread_join(threads[i], NULL);
		}

		free(minb);
		fprintf(stderr, "KDE: h = %lf, nSet = %d\n", h, nSet);
		fflush(stderr);
		for (int i = 0; i < (int)((nbStdDev * 2 * (maxi_h + 1) + 1)); i++)
			free(gaussian[i]);
		free(gaussian);
	}

	return 0;
};

int KDE::applyGaussianFilter()
{
	int startRow, startCol, endRow, endCol;
	while (true)
	{
		if (job.empty())
		{
			return 0;
		}
		pthread_mutex_lock(&kde_mutex);
		if (job.empty())
		{
			pthread_mutex_unlock(&kde_mutex);
			return 0;
		}
		startRow = job.front().first;
		startCol = job.front().second;
		job.pop();

		pthread_mutex_unlock(&kde_mutex);
		endRow = min(O.endRow(), startRow + this->rowStep);
		endCol = startCol + this->colStep;
		for (int i = startRow; i < endRow; i++)
		{
			for (int j = max(startCol, O.startCol(i));
				 j < min(endCol, O.endCol(i)); j++)
			{
				if (T.get(i, j) < 0)
				{
					// check if there are enough entries in the square around i,j
					int lowi =
						(i - nbStdDevTh) >= O.startRow() ? (i - nbStdDevTh) : O.startRow();
					int highi =
						(i + nbStdDevTh) < O.endRow() ? (i + nbStdDevTh) : O.endRow() - 1;

					int lowj =
						(j - nbStdDevTh) >= O.startCol(lowi) ? (j - nbStdDevTh) : O.startCol(lowi);
					int highj =
						(j + nbStdDevTh) < O.endCol(lowi) ? (j + nbStdDevTh) : O.endCol(lowi) - 1;

					int nHits = 0;
					if (lowi > O.startRow() && lowj > O.startCol(O.startRow()))
					{
						nHits = cumO->get(highi, highj) - cumO->get(lowi - 1, highj);
						if (highi <= lowj - 1)
						{
							nHits -= cumO->get(highi, lowj - 1);
						}
						else
						{
							nHits -= cumO->get(lowj - 1, lowj - 1);
						}

						if (lowi - 1 <= lowj - 1)
						{
							nHits += cumO->get(lowi - 1, lowj - 1);
						}
						else
						{
							nHits += cumO->get(lowj - 1, lowj - 1);
						}
					}

					if (lowi > O.startRow() && lowj == O.startCol(O.startRow()))
						nHits = cumO->get(highi, highj) - cumO->get(lowi - 1, highj);
					if (lowi == O.startRow() && lowj > O.startCol(O.startRow()))
					{
						nHits = cumO->get(highi, highj);
						if (highi <= lowj - 1)
							nHits -= cumO->get(highi, lowj - 1);
						else
							nHits -= cumO->get(lowj - 1, lowj - 1);
					}
					if (lowi == O.startRow() && lowj == O.startCol(O.startRow()))
						nHits = cumO->get(highi, highj);

					if ((lastIter && nHits > 0) || nHits >= option_kdeMinCount)
					{
						//	      fprintf(stderr,"computing %d %d nHits=%d\n",i,j,nHits);

						double s = 0;
						int c = 0;
						double sumw = 0;
						double gs = gaussianSum;
						bool isSafe;
						// check if entire square around i,j is away from diagonal and borders
						if (i - maxRadius >= O.startRow() && i + maxRadius < O.endRow() && j + maxRadius < O.endCol(i + maxRadius) && j - maxRadius >= O.startCol(i - maxRadius) && i + maxRadius < j - maxRadius)
							isSafe = true;
						else
							isSafe = false;

						if (isSafe)
						{
							int mr = (int)maxRadius;
							for (int ii = i - mr; ii <= i + mr; ii++)
							{
								int range = minb[ii - i + mr];
								int jj = j - range;
								if (O.get(ii, jj) == 0)
									jj = nextNonEmpty->get(ii, jj);

								while (jj != -1 && jj <= j + mr)
								{
									if (ii == jj)
										continue; // skip items on main diagonal
									if (jj - j <= range)
									{
										s +=
											OoverBias.get(ii, jj) * gaussian[ii - i + nbStdDevTh][jj - j + nbStdDevTh];
									}
									jj = nextNonEmpty->get(ii, jj);
								}
							}
						}
						else
						{
							for (int a = -maxRadius; a <= maxRadius; a++)
							{
								int ii = i + a;
								int range = minb[(int)(a + maxRadius)];
								for (int b = -range; b <= range; b++)
								{
									int jj = j + b;
									if (ii != jj && (ii >= O.startRow() && ii < O.endRow() && jj >= ii && jj >= O.startCol(ii) && jj < O.endCol(ii)))
									{
										float x = OoverBias.get(ii, jj);
										if (x > 0)
										{
											s += x * gaussian[a + nbStdDevTh][b + nbStdDevTh];
											//			fprintf(stderr,"Adding_s %d %d %f %f\n",ii,jj,OoverBias.get(ii,jj),gaussian[a+nbStdDevTh][b+nbStdDevTh]);
										}
									}
									else
									{
										gs -= gaussian[a + nbStdDevTh][b + nbStdDevTh];
									}
								}
							}
						}
						//	      fprintf(stderr,"Setting T[%d][%d] = %f\n",i,j,s/gs);
						//						if(i!=j)
						T.set(i, j, s / gs);

						nSet++;
					}
				}
				else
				{
					nSet++;
				}
			}
		}
	}
	return 0;
}
