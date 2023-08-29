/*
 * MRF.cpp
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#include "MRF.h"
#include "utils.h"
#include "readOptions.h"
#include "optimize.h"

pthread_mutex_t mrf_mutex;

int MRF::initialize()
{
	rowStep = 500;
	colStep = 500;
	setBoundaries();
	fprintf(stderr, "Done setting up boundaries\n");

	hasChanged = new UpperDiag<char>(option_firstRow, option_firstCol,
									 option_lastRow, option_lastCol);

	// initalization of T
	if (!(option_kdeInitialization || method_kde))
	{
		for (int i = O.startRow(); i < O.endRow(); i++)
		{
			for (int j = O.startCol(i); j < O.endCol(i); j++)
			{
				T.set(i, j, O.get(i, j) / bias[i] / bias[j] + 0.001);
				// T.set(i, j, drand48());
				cout << "setting\n";
			}
		}
	}

	double originalLike = evaluateFullLikelihood();
	cout << "MRF Original Likelihood: " << originalLike << endl;

	nChanges = 1;

	return 0;
}
int MRF::runMRF()
{
	int tids[n_threads];
	pthread_t threads[n_threads];
	pthread_mutex_init(&mrf_mutex, NULL);
	for (rep = 0; rep < option_mrfMaxIter; rep++)
	{
		if (option_testMatrix[0])
			evaluateAgainstTest();

		cout << "MRF: " << rep << endl;

		double sumChanges = 0;

		for (int i = O.startRow(); i < O.endRow(); i += rowStep)
		{
			for (int j = O.startCol(i); j < O.endCol(i); j += colStep)
			{
				job.push(make_pair(i, j));
			}
		}

		for (int i = 0; i < n_threads; i++)
		{
			tids[i] = i;
			pthread_create(&threads[i], NULL, blockMRF, this);
		}

		for (int i = 0; i < n_threads; i++)
		{
			pthread_join(threads[i], NULL);
		}
		nChanges = 0;
		sumChanges = 0;

		fprintf(stderr, "rep=%d: nChanges=%ld, sumCHanges=%lf\n", rep, nChanges,
				sumChanges);

		fflush(stderr);
		cout << "MRF iter: " << rep + 1 << " Likelihood: "
			 << evaluateFullLikelihood() << endl
			 << flush;

	} // end of MRF iteration loop
}

int MRF::blockMRF()
{
	int startRow, startCol, endRow, endCol;
	double *allNei, newT;
	allNei = new double[9];
	int nNei;
	double t, o;
	while (true)
	{
		if (job.empty())
		{
			return 0;
		}
		pthread_mutex_lock(&mrf_mutex);
		if (job.empty())
		{
			pthread_mutex_unlock(&mrf_mutex);
			return 0;
		}
		startRow = job.front().first;
		startCol = job.front().second;
		job.pop();
		pthread_mutex_unlock(&mrf_mutex);
		endRow = min(O.endRow(), startRow + this->rowStep);
		endCol = startCol + this->colStep;

		int neiLeft, neiRight, neiBottom, neiTop;
		bool changed;
		for (int i = startRow; i < endRow; i++)
		{
			for (int j = max(startCol, O.startCol(i));
				 j < min(endCol, O.endCol(i)); j++)
			{
				if (i == j) //skip self-interaction
					continue;
				t = T.get(i, j);
				o = O.get(i, j);
				getNeighbors(i, j, nNei, allNei);

				newT = greedyOptimizeT(o, t, allNei, nNei, bias[i] * bias[j]);
				T.set(i, j, newT);
			}
		} // end of for i for j
	}
	return 0;
}
