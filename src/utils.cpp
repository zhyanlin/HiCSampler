/*
 * utils.cpp
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#include "utils.h"
#include "readOptions.h"
#include <fstream>
#include "opt_median.h"

 annotation anno;

float factor[NFACTORS] = {0.5, 0.7, 0.9, 0.95, 0.97, 0.99, 1.0, 1.01, 1.03,
						  1.05, 1.1, 1.3, 1.5};
float logFactor[NFACTORS];

int fullMatrix_firstRow, fullMatrix_firstCol, fullMatrix_lastRow,
	fullMatrix_lastCol, fullMatrix_size;

int nNeighborhoods = 11;

int chr1, chr2;

bool fast = true;

float *bias, *piror_k, *piror_theta;

int numberOfDumpables;

float *logKFact;

double *precomputedNormPDF;
double *precomputedNormLogPDF;

int fullRange = 100 * NORM_RES;
int floatFullRange = 2 * 100 * NORM_RES;

void outputSparseBoundaryMatrix(char *fn);

void setBoundaries()
{

	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i); j < O.endCol(i); j++)
		{
			int o = O.get(i, j);
			if (o > 0)
			{
				float x = o / (bias[i] * bias[j]);
				OoverBias.set(i, j, x);
			}
		}
	}

	float sqrtLeftOver2[100000];
	for (int i = 0; i < 100000; i++)
		sqrtLeftOver2[i] = sqrt(i / 2.0);

	// horizontal
	for (int i = O.startRow(); i < O.endRow() - 1; i++)
	{

		// consider a boundary starting down from (i,j)
		int distLeft[10000];
		int distRight[10000];
		int maxj = -1;
		int maxj2 = -1;
		float maxks = -1;

		for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
		{

			if (OoverBias.get(i, j) == 0 && j + 1 >= O.startCol(i) && j + 1 < O.endCol(i) && OoverBias.get(i + 1, j) == 0)
				continue;

			memset(distLeft, 0, 10000 * sizeof(int));
			memset(distRight, 0, 10000 * sizeof(int));
			int nLeft = 0;
			int nRight = 0;
			int nbNonMax = 0;
			float maxSearch = 0;
			int maxEntry = 0;

			int j2 = 0;

			for (j2 = j; j2 < O.endCol(i) && nbNonMax < MAX_NBNONMAX; j2++)
			{ // end point of boundary
				nbNonMax++;
				// only consider end points where there is data
				if (horizontalBoundary.get(i, j2))
					break;
				int entryLeft = (int)(OoverBias.get(i, j2) + 0.49);
				int entryRight = (int)(OoverBias.get(i + 1, j2) + 0.49);

				distLeft[entryLeft]++;
				nLeft++;
				distRight[entryRight]++;
				nRight++;
				if (entryLeft > maxEntry)
					maxEntry = entryLeft;
				if (entryRight > maxEntry)
					maxEntry = entryRight;

				if (entryLeft == 0 && entryRight == 0)
					continue;

				int cumLeft2 = 0;
				int cumRight2 = 0;
				int maxDiff2 = 0;
				for (int a = 0; a < maxEntry; a++)
				{
					cumLeft2 += distLeft[a];
					cumRight2 += distRight[a];
					if (abs(cumLeft2 - cumRight2) > maxDiff2)
					{
						maxDiff2 = abs(cumLeft2 - cumRight2);
					}
				}

				float ks = ((float)maxDiff2 / nLeft) * sqrtLeftOver2[nLeft];

				if (ks > maxSearch)
				{
					nbNonMax = 0;
					maxSearch = ks;
				}

				if (ks > maxks)
				{

					maxks = ks;
					maxj = j;
					maxj2 = j2;
				}
				else
				{
					nbNonMax++;
				}
			}
		}

		if (maxks > option_boundaryKS)
		{
			// fprintf(stderr, "For row %d, boundary from %d to %d, with ks=%lf\n",
			// 		i, maxj, maxj2, maxks);
			// fflush(stderr);
			for (int j = maxj; j <= maxj2; j++)
			{
				horizontalBoundary.set(i, j, 1);
			}
			i--;
		}
	}

	// vertical
	for (int j = O.startCol(O.startRow()); j < O.endCol(O.startRow()) - 1;
		 j++)
	{

		// consider a boundary starting down from (i,j)
		int distLeft[10000];
		int distRight[10000];
		int maxi = -1;
		int maxi2 = -1;
		float maxks = -1;

		for (int i = O.startRow(); i < O.endRow(); i++)
		{
			if (i >= j)
				continue;
			if (OoverBias.get(i, j) == 0 && OoverBias.get(i, j + 1) == 0)
				continue;

			memset(distLeft, 0, 10000 * sizeof(int));
			memset(distRight, 0, 10000 * sizeof(int));
			int nLeft = 0;
			int nRight = 0;
			int nbNonMax = 0;
			float maxSearch = 0;
			int maxEntry = 0;

			int i2 = 0;

			for (i2 = i; i2 <= j && i2 < O.endRow() && nbNonMax < MAX_NBNONMAX;
				 i2++)
			{ // end point of boundary
				nbNonMax++;
				// only consider end points where there is data
				if (verticalBoundary.get(i2, j))
					break;
				int entryLeft = (int)(OoverBias.get(i2, j) + 0.49);
				int entryRight = (int)(OoverBias.get(i2, j + 1) + 0.49);

				distLeft[entryLeft]++;
				nLeft++;
				distRight[entryRight]++;
				nRight++;
				if (entryLeft > maxEntry)
					maxEntry = entryLeft;
				if (entryRight > maxEntry)
					maxEntry = entryRight;

				if (entryLeft == 0 && entryRight == 0)
					continue;

				int cumLeft2 = 0;
				int cumRight2 = 0;
				int maxDiff2 = 0;
				for (int a = 0; a < maxEntry; a++)
				{
					cumLeft2 += distLeft[a];
					cumRight2 += distRight[a];
					if (abs(cumLeft2 - cumRight2) > maxDiff2)
					{
						maxDiff2 = abs(cumLeft2 - cumRight2);
					}
				}

				float ks = ((float)maxDiff2 / nLeft) * sqrtLeftOver2[nLeft];

				if (ks > maxSearch)
				{
					nbNonMax = 0;
					maxSearch = ks;
				}

				if (ks > maxks)
				{

					maxks = ks;
					maxi = i;
					maxi2 = i2;
				}
				else
				{
					nbNonMax++;
				}
			}
		}

		if (maxks > option_boundaryKS)
		{
			// fprintf(stderr,
			// 		"For columns %d, vertical boundary from %d to %d, with ks=%lf\n",
			// 		j, maxi, maxi2, maxks);
			// fflush(stderr);
			for (int i = maxi; i <= maxi2; i++)
			{
				verticalBoundary.set(i, j, 1);
			}
			j--;
		}
	}

	// additional boundaries for large un-called regions
	int *sumOfRow, *sumOfCol;
	sumOfRow = new int[O.endRow()]();
	sumOfCol = new int[O.endCol(O.endRow())]();
	memset(sumOfRow, 0, O.endRow() * sizeof(int));
	memset(sumOfCol, 0, O.endCol(O.endRow()) * sizeof(int));
	for (int i = O.startRow(); i < O.endRow() - 1; i++)
	{
		for (int j = O.startCol(i); j < O.endCol(i); j++)
		{
			sumOfRow[i] += O.get(i, j);
			sumOfCol[j] += O.get(i, j);
		}
	}
	bool isLargeUnCalledRegionBoundary = false;
	int uncalledRF = 0;
	//horizontal
	for (int i = O.startRow(); i < O.endRow() - 1; i++)
	{
		if (sumOfRow[i] > 0)
		{
			if (uncalledRF > 4 && i > 0)
			{
				for (int j = O.startCol(i - 1); j < O.endCol(i - 1); j++)
					horizontalBoundary.set(i - 1, j, 1);
			}
			uncalledRF = 0;
			isLargeUnCalledRegionBoundary = true;
			for (int ii = 1; ii <= 5 && i + ii < O.endRow(); ii++)
			{
				if (sumOfRow[i + ii] > 0)
				{
					isLargeUnCalledRegionBoundary = false;
					continue;
				}
			}
			if (isLargeUnCalledRegionBoundary)
			{
				for (int j = O.startCol(i); j < O.endCol(i); j++)
					horizontalBoundary.set(i, j, 1);
			}
		}
		else
			uncalledRF += 1;
	}
	//vertical
	for (int j = O.startCol(O.startRow()); j < O.endCol(O.endRow()); j++)
	{
		if (sumOfCol[j] > 0)
		{
			if (uncalledRF > 4 && j > 0)
			{
				for (int i = O.startRow(); i < O.endRow(); i++)
				{
					if (j - 1 > O.startCol(i) && j - 1 < O.endCol(i))
						verticalBoundary.set(i, j - 1, 1);
				}
			}
			uncalledRF = 0;
			isLargeUnCalledRegionBoundary = true;
			for (int jj = 1; jj <= 5 && j + jj < O.endCol(O.endRow()); jj++)
			{
				if (sumOfCol[j + jj] > 0)
				{
					isLargeUnCalledRegionBoundary = false;
					continue;
				}
			}
			if (isLargeUnCalledRegionBoundary)
			{
				for (int i = O.startRow(); i < O.endRow(); i++)
				{
					if (j > O.startCol(i) && j < O.endCol(i))
						verticalBoundary.set(i, j, 1);
				}
			}
		}
		else
			uncalledRF += 1;
	}

	// end of additional boundaries for large un-called regions

	// thinning verticalBoundary
	int tmp;
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i); j < O.endCol(i) - 1; j++)
		{
			if ((int)verticalBoundary.get(i, j) > 0)
			{
				if ((int)verticalBoundary.get(i, j + 1) > 0)
				{
					verticalBoundary.set(i, j + 1, verticalBoundary.get(i, j + 1) + verticalBoundary.get(i, j));
					verticalBoundary.set(i, j, 0);
				}
				else
				{
					tmp = (int)verticalBoundary.get(i, j);
					tmp = tmp / 2;
					verticalBoundary.set(i, j, 0);
					verticalBoundary.set(i, j - tmp, 1);
					// cout << "I'm a boundary " << i << " " << j - tmp << " " << 1 << endl;
				}
			}
		}
	}

	// thinning horizontalBoundary
	for (int i = O.startRow(); i < O.endRow() - 1; i++)
	{
		for (int j = O.startCol(i); j < O.endCol(i); j++)
		{
			if ((int)horizontalBoundary.get(i, j) > 0)
			{
				if ((int)horizontalBoundary.get(i + 1, j) > 0)
				{
					horizontalBoundary.set(i + 1, j, horizontalBoundary.get(i + 1, j) + horizontalBoundary.get(i, j));
					horizontalBoundary.set(i, j, 0);
				}
				else
				{
					tmp = (int)horizontalBoundary.get(i, j);
					tmp = tmp / 2;
					horizontalBoundary.set(i, j, 0);
					horizontalBoundary.set(i - tmp, j, 1);
					// cout << "I'm a boundary " << i - tmp << " " << j << " " << 1 << endl;
				}
			}
		}
	}

	if (option_boundaryOutput[0])
	{
		outputSparseBoundaryMatrix(option_boundaryOutput);
	}
}

// returns the median (or mean) of the neighborhood of i,j
float getMedianNeighbor(int i, int j)
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

				allNei[nNei] = T.get(ii, jj); // + drand48() * 0.0001;
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
// flag for debugging
int toPrint;

double evaluateAgainstTest()
{
	// get sum of T and test, to normalize
	double sumO = 0;
	double sumTestO = 0;
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
		{
			sumO += O.get(i, j); // * bias[i] * bias[j];
			sumTestO += testO.get(i, j);
		}
	}

	float scale = sumTestO / sumO;
	double SSE = 0;
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
		{
			SSE += pow(
				testO.get(i, j) - T.get(i, j) * bias[i] * bias[j] * scale,
				2);
		}
	}
	fprintf(stderr, "SSE = %lf\n", SSE);
	return SSE;
}

double evaluateAgainstTest(UpperDiag<float> *currentT)
{
	// get sum of T and test, to normalize
	double sumO = 0;
	double sumTestO = 0;
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
		{
			sumO += O.get(i, j); // * bias[i] * bias[j];
			sumTestO += testO.get(i, j);
		}
	}

	float scale = sumTestO / sumO;

	double SSE = 0;
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
		{

			SSE += pow(
				testO.get(i, j) - currentT->get(i, j) * bias[i] * bias[j] * scale,
				2);
		}
	}

	return SSE;
}

void precomputeStuff()
{
	// normal distribution PDF
	precomputedNormPDF = (double *)malloc((2000000 + 1) * sizeof(double));
	precomputedNormLogPDF = (double *)malloc((2000000 + 1) * sizeof(double));
	double s = 0;
	float denom = sqrt(6.28318530717959);
	for (int i = -NORM_RES * 100; i <= NORM_RES * 100; i++)
	{

		precomputedNormPDF[i + NORM_RES * 100] = exp(
													 -((float)i / NORM_RES * (float)i / NORM_RES) / 2) /
												 denom;
		precomputedNormLogPDF[i + NORM_RES * 100] = -((float)i / NORM_RES * (float)i / NORM_RES) / 2 - log(denom);

		s += precomputedNormPDF[i + NORM_RES * 100];
	}

	for (int i = -NORM_RES * 100; i <= NORM_RES * 100; i++)
	{
		precomputedNormPDF[i + NORM_RES * 100] /= s;
		precomputedNormLogPDF[i + NORM_RES * 100] -= log(s);
	}

	// log of multiplicative factors
	for (int i = 0; i < NFACTORS; i++)
		logFactor[i] = log(factor[i]);
}

void allocateMemory()
{
	fprintf(stderr, "Alloc %d %d %d %d", option_firstRow, option_firstCol,
			option_lastRow, option_lastCol);

	anno.gc=new double[option_lastRow-option_firstRow+1];
	anno.len=new double[option_lastRow-option_firstRow+1];
	anno.map=new double[option_lastRow-option_firstRow+1];
	for(int i=0;i<option_lastRow-option_firstRow+1;i++){
		anno.gc[i]=0;
		anno.len[i]=0;
		anno.map[i]=0;
	}
	T = UpperDiag<float>(option_firstRow, option_firstCol, option_lastRow,
						 option_lastCol);
	O = UpperDiag<float>(option_firstRow, option_firstCol, option_lastRow,
						 option_lastCol);
	if (option_testMatrix[0])
		testO = UpperDiag<float>(option_firstRow, option_firstCol,
								 option_lastRow, option_lastCol);
	OoverBias = UpperDiag<float>(option_firstRow, option_firstCol,
								 option_lastRow, option_lastCol);

	verticalBoundary = UpperDiag<char>(option_firstRow, option_firstCol,
									   option_lastRow, option_lastCol);
	horizontalBoundary = UpperDiag<char>(option_firstRow, option_firstCol,
										 option_lastRow, option_lastCol);

	bias = (float *)malloc(fullMatrix_size * sizeof(float));
	piror_k = (float *)malloc(60000 * sizeof(float));
	piror_theta = (float *)malloc(60000 * sizeof(float));
	for (int i = 0; i < fullMatrix_size; i++)
	{
		bias[i] = 1;
		piror_k[i] = 0.01;
		piror_theta[i] = 0.01;
	}
}

void getMatrixSize(char *fn)
{

	fullMatrix_firstRow = fullMatrix_firstCol = 999999999;
	fullMatrix_lastRow = fullMatrix_lastCol = -1;
	FILE *f = fopen(fn, "r");

	if (option_inputSparseMatrix)
	{
		char line[100000];
		fprintf(stderr, "Reading input sparse input file to get sizes\n");

		char c1[100];
		char c2[100];

		while (fscanf(f, "%[^\n]\n", line) != EOF)
		{
			// if (line[0] == '#'){
			
			
			// sscanf(line,"#	%d	%d	%d	%d", &fullMatrix_firstRow,&fullMatrix_lastRow,&fullMatrix_firstCol,&fullMatrix_lastCol);
			
			// 	continue;
				
				
			// }
			float foo;
			int p1, p2;
			sscanf(line, "%s %d %s %d	%f", c1, &p1, c2, &p2,&foo);
			if (p1 < fullMatrix_firstRow)
				fullMatrix_firstRow = p1;
			if (p2 < fullMatrix_firstCol)
				fullMatrix_firstCol = p2;

			if (p1 > fullMatrix_lastRow)
				fullMatrix_lastRow = p1;
			if (p2 > fullMatrix_lastCol)
				fullMatrix_lastCol = p2;
		}
		fullMatrix_lastCol+=1;
		fullMatrix_lastRow+=1;

		if (fullMatrix_firstCol > fullMatrix_firstRow)
			fullMatrix_firstCol = fullMatrix_firstRow;
		if (fullMatrix_lastCol > fullMatrix_lastRow)
			fullMatrix_lastRow = fullMatrix_lastCol;

		fprintf(stderr, "firstRow=%d, firstCol=%d, lastRow=%d, lastCol=%d\n",
				fullMatrix_firstRow, fullMatrix_firstCol, fullMatrix_lastRow,
				fullMatrix_lastCol);
		fflush(stderr);
	}

	// fullMatrix_size = max(fullMatrix_lastRow, fullMatrix_lastCol) + 1;
	fullMatrix_size = fullMatrix_lastRow + 1; // 0-start

	if (option_firstRow == -1)
	{
		option_firstRow = fullMatrix_firstRow; //0;
		option_firstCol = fullMatrix_firstCol; //0;
		option_lastRow = fullMatrix_lastRow;
		option_lastCol = fullMatrix_lastCol;
	}

	fclose(f);

	fflush(stderr);
}

void readFullOMatrix(char *fn, UpperDiag<float> *mat, bool setBias)
{

	FILE *f = fopen(fn, "r");
	char line[1000];

	double fullSum = 0;
	double *sumRow = (double *)malloc(fullMatrix_size * sizeof(double));
	for (int i = 0; i < fullMatrix_size; i++)
	{
		sumRow[i] = BIAS_PSEUDOCOUNT;
	}
	fullSum += fullMatrix_size * BIAS_PSEUDOCOUNT;

	if (option_inputSparseMatrix)
	{
		while (fscanf(f, "%[^\n]\n", line) != EOF)
		{
			if (line[0] == '#')
				continue;

			int foo;
			int p1, p2;
			float c;
			char c1[100];
			char c2[100];
			sscanf(line, "%s %d %s %d	%f",  c1, &p1, c2, &p2,&c);

			sumRow[p1] += c;
			sumRow[p2] += c;
			fullSum += c + c;
			if(abs(p2-p1)<bandSize)
			if (!option_maxFragDistance || (abs(p2 - p1) <= option_maxFragDistance))
			{
				if (option_lastRow == -1 || (p1 >= option_firstRow && p1 <= option_lastRow && p2 >= option_firstCol && p2 <= option_lastCol))
				{
					if (p1 < option_firstRow || p1 > option_lastRow)
						continue;
					if (p2 < option_firstCol || p1 > option_lastCol)
						continue;
					mat->set(p1, p2, c);
				}
			}
		}
	}
	fclose(f);

	fprintf(stderr, "fullSum=%f\n", fullSum);

	if (setBias)
	{
		for (int i = 0; i < fullMatrix_size; i++)
		{
			bias[i] = sumRow[i] / ((fullSum) / fullMatrix_size);
			if (bias[i] > MAX_BIAS)
				bias[i] = MAX_BIAS;
			if (bias[i] < 1.0 / MAX_BIAS)
				bias[i] = 1.0 / MAX_BIAS;
			if (bias[i] < 1)
				bias[i] = 1.0;
			// fprintf(stderr, "Bias[%d]=%f\n", i, bias[i]);
		}

		if (option_bias)
		{
			for (int i = 0; i < fullMatrix_size; i++)
				bias[i] = 0;
			std::ifstream biasFile(option_bias);
			int pos;
			float bv, minBias;
			minBias = 0.1;
			while (biasFile >> pos >> bv)
			{
				if (pos < fullMatrix_size)
				{
					bias[pos] = bv;
					fprintf(stderr, "Bias[%d]=%f\n", pos, bias[pos]);
				}
			}
			for (int i = 0; i < fullMatrix_size; i++)
				if (bias[i] < minBias)
					bias[i] = 1;
		}
	}
	//read piror
	if(option_prior_file){
	std::ifstream pirorFile(option_prior_file);
	int genomicDistance;
	float pirork, pirortheta;
	while (pirorFile >> genomicDistance >> pirork >> pirortheta)
	{
		piror_k[genomicDistance] = pirork;
		piror_theta[genomicDistance] = pirortheta;
	}
	cout << "finished loading piror\n"
		 << flush;
		 }
	//end of read piror
}

void computeTMatrix_fixed()
{
	cout << "computeTMatrix_fixed in progress...\n";
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i); j < O.endCol(i); j++)
		{
			T.set(i, j, O.get(i, j) / (bias[i] * bias[j]));	
		}
	}
	cout<<"end of computeTMatrix_fixed\n";
}

void computeBias()
{
	if (option_noBias)
	{
		for (int i = 0; i < fullMatrix_size; i++)
		{
			bias[i] = 1;
		}
		return;
	}

	// calculate biases
	fprintf(stderr, "ComputeBias\n");
	fflush(stderr);
	double fullSum = 0;
	double *sumRow = (double *)malloc(fullMatrix_size * sizeof(double));

	for (int i = 0; i < fullMatrix_size; i++)
	{
		sumRow[i] = BIAS_PSEUDOCOUNT;
	}
	fullSum += fullMatrix_size * BIAS_PSEUDOCOUNT;

	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i); j < O.endCol(i); j++)
		{
			int x = O.get(i, j);
			sumRow[i] += x;
			sumRow[j] += x;
			fullSum += x + x;
		}
	}

	fprintf(stderr, "fullSum=%f\n", fullSum);

	for (int i = 0; i < fullMatrix_size; i++)
	{
		bias[i] = sumRow[i] / ((fullSum) / fullMatrix_size);
		if (bias[i] > MAX_BIAS)
			bias[i] = MAX_BIAS;
		if (bias[i] < 1.0 / MAX_BIAS)
			bias[i] = 1.0 / MAX_BIAS;
		fprintf(stderr, "Bias[%d]=%f\n", i, bias[i]);
	}

	free(sumRow);
}

void outputSparseMatrix(char *fn)
{
	FILE *out = fopen(fn, "w");
	if (option_lastRow != -1 && option_lastCol != -1)
		fprintf(out, "# %d %d %d %d\n", option_firstRow, option_lastRow,
				option_firstCol, option_lastCol);
	else
		fprintf(out, "# %d %d %d %d\n", option_firstRow, option_lastRow,
				option_firstCol, option_lastCol);

	for (int i = T.startRow(); i < T.endRow(); i++)
	{

		for (int j = T.startCol(i); j < T.endCol(i); j++)
		{
			float x = T.get(i, j);
			if (!option_outputNormalized)
				x *= (bias[i] * bias[j]);
			if (x > option_minOutput)
				fprintf(out, "%5.3lf %d %d %d %d\n", x / option_minOutput, chr1,
						i, chr2, j);
		}
	}
	fclose(out);
}

void outputSparseMatrixGZ(char *fn)
{
	gzFile gzO = gzopen(fn, "w");

	if (option_lastRow != -1 && option_lastCol != -1)
		gzprintf(gzO, "# %d %d %d %d\n", option_firstRow, option_lastRow,
				 option_firstCol, option_lastCol);
	else
		gzprintf(gzO, "# %d %d %d %d\n", option_firstRow, option_lastRow,
				 option_firstCol, option_lastCol);

	for (int i = T.startRow(); i < T.endRow(); i++)
	{

		for (int j = T.startCol(i); j < T.endCol(i); j++)
		{
			if (i == j)
				continue;
			float x = T.get(i, j);
			// if (!option_outputNormalized)
			// 	x *= (bias[i] * bias[j]);
			// if (x > option_minOutput)
			gzprintf(gzO, "%d %d %d %d %5.3lf\n",
					 chr1, i, chr2, j, x);
		}
	}
	gzclose(gzO);
}

void outputSparseBoundaryMatrix(char *fn)
{
	cout << fn << endl;
	// return;
	FILE *out = fopen(fn, "w");
	// return;
	if (option_lastRow != -1 && option_lastCol != -1)
		fprintf(out, "# %d %d %d %d\n", option_firstRow, option_lastRow,
				option_firstCol, option_lastCol);
	else
		fprintf(out, "# %d %d %d %d\n", option_firstRow, option_lastRow,
				option_firstCol, option_lastCol);
	// return ;
	for (int i = horizontalBoundary.startRow(); i < horizontalBoundary.endRow();
		 i++)
	{
		for (int j = horizontalBoundary.startCol(i) + 1;
			 j < horizontalBoundary.endCol(i) - 1; j++)
		{
			// horizontalBoundary.get(i, j);
			// 	verticalBoundary.get(i, j);
			// cout<<i<<","<<j<<endl;
			if (horizontalBoundary.get(i, j)
				// || (i != horizontalBoundary.startRow() && horizontalBoundary.get(i - 1, j))
				|| verticalBoundary.get(i, j)
				// || (j != horizontalBoundary.startCol(i) + 1 && verticalBoundary.get(i, j - 1))
			)
				fprintf(out, "1 %d %d %d %d\n", chr1, i, chr2,
						j);
		}
	}
	fclose(out);
}

double norm_Randn(float mu, float sigma)
{
	double epsilon = 1e-10;
	double two_pi = 2.0 * 3.14159265358979323846;
	double z0;
	double u1, u2;
	do
	{
		u1 = drand48();
		u2 = drand48();
	} while (u1 <= epsilon);
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	return z0 * sigma + mu;
};

float normpdf(const float &x, const float &m, const float &s)
{
	int i = ((x - m) / s) * NORM_RES + fullRange;
	if (i < 0 || i > floatFullRange)
		return 0;
	return precomputedNormPDF[i];
};

void precomputeLogKFact()
{
	logKFact = (float *)malloc(1000000 * sizeof(float));

	logKFact[0] = 0;
	for (int i = 1; i < 1000000; i++)
	{
		logKFact[i] = logKFact[i - 1] + log(i);
	}
};

double median(int nNei, float *allNei)
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
};
