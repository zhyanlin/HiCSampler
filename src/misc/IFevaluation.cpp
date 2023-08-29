#include "../utils.h"
#include "../contactMap.h"
using namespace std;

double evaluateFullLikelihood()
{
	double like = 0, l;
	double lModel = 0;
	double lObsGivenModel = 0;
	float t, medianNeighbor;
	int o;
	double biases;
	double variance = 1;
	double ld, lm;
	for (int i = O.startRow(); i < O.endRow(); i++)
	{
		for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
		{
			ld = 0;
			lm = 0;
			t = T.get(i, j);
			o = O.get(i, j);
			biases = bias[i] * bias[j];
			 ld= -biases * t - logKFact[o] + o * log(biases * t + 0.00001);
			 lObsGivenModel+=ld;

			if (j + 1 < O.endCol(i) && !verticalBoundary.get(i, j))
				lm -= pow(log(1 + T.get(i, j + 1)) - log(1 + t), 2) / variance;
			if (i + 1 < O.endRow())
			{
				if (j > O.startCol(i) + 1 && !horizontalBoundary.get(i, j))
					lm -= pow(log(1 + T.get(i + 1, j)) - log(1 + t), 2) / variance;
				if (j + 1 < O.endCol(i + 1))
					if ((!horizontalBoundary.get(i, j) && !verticalBoundary.get(i + 1, j)) ||
						(!verticalBoundary.get(i, j) && !horizontalBoundary.get(i, j + 1)))
						lm -= pow(log(1 + T.get(i + 1, j + 1)) - log(1 + t), 2) / variance;
				if (j - 1 > O.startCol(i + 1))
					if ((!verticalBoundary.get(i, j - 1) && !horizontalBoundary.get(i, j - 1)) ||
						(!horizontalBoundary.get(i, j) && !verticalBoundary.get(i + 1, j - 1)))
						lm -= pow(log(1 + T.get(i + 1, j - 1)) - log(1 + t), 2) / variance;
			}
			lModel+=lm;
			cout<<i<<","<<j<<": "<<o<<" "<<t<<" "<<biases<<" "<<ld<<" "<<lm<<endl;
		}
	}
	cout<<"full log-likelihood= "<< lObsGivenModel + lModel <<endl;
	return lObsGivenModel + lModel;
}

// double evaluateFullLikelihoodSkipBoundary()
// {
// 	double like = 0;

// 	for (int i = O.startRow(); i < O.endRow(); i++)
// 	{
// 		for (int j = O.startCol(i) + 1; j < O.endCol(i); j++)
// 		{
// 			bool skip = false;
// 			int ii, jj, off_j;
// 			int startRow = O.startRow(), endRow = O.endRow();
// 			for (int off_i = -1; off_i <= 1; off_i++)
// 			{
// 				for (off_j = -1; off_j <= 1; off_j++)
// 				{
// 					ii = i + off_i;
// 					jj = j + off_j;
// 					if (ii < jj && ii >= startRow && ii < endRow && jj >= O.startCol(ii) && jj < O.endCol(ii))
// 					{

// 						if (off_j == -1 && verticalBoundary.get(i, j - 1))
// 							skip = true;
// 						if (off_j == +1 && verticalBoundary.get(i, j))
// 							skip = true;
// 						if (off_i == -1 && (i == startRow || horizontalBoundary.get(i - 1, j)))
// 							skip = true;
// 						if (off_i == +1 && horizontalBoundary.get(i, j))
// 							skip = true;
// 					}
// 				}
// 			}
// 			if (skip)
// 				continue;
// 			float medianNeighbor = getMedianNeighbor(i, j);
// 			float t = T.get(i, j);
// 			int o = O.get(i, j);
// 			double l = evaluateLikelihood(t, o, medianNeighbor, bias[i] * bias[j]);
// 			like += l;
// 			cout << i << " " << j << " " << l << endl;
// 		}
// 	}
// 	cout << "FullLikelihood= " << like << endl;
// 	return like;
// }
void readFullTMatrix(char *fn, UpperDiag<float> *mat)
{
	gzFile gzIn = gzopen(fn, "r");
	char line[1000];
	int foo;
	int p1, p2;
	float c;
	char c1[100];
	char c2[100];

	for (int i = T.startRow(); i < T.endRow(); i++)
	{

		for (int j = T.startCol(i); j < T.endCol(i); j++)
		{
			mat->set(i, j, 0);
		}
	}

	while (!gzeof(gzIn))
	{
		if (gzgets(gzIn, line, 1000))
		{
			sscanf(line, "%s %d %s %d %f", c1, &p1, c2, &p2, &c);
			mat->set(p1, p2, c);
		}
	}
};
int main(int argc, char *argv[])
{

	srand48(time(NULL));
	if (argc != 3 && argc != 4)
	{
		fprintf(stderr,
				"IFevaluation <matrix_file_with_row_count> <matrix_file_with_IF>  0/1[skipBoundary optional,default is 0.this is not implemented]\n");
		exit(1);
	}

	if (option_gibbs)
	{ // gibbs is not implemented yet
		exit(1);
	}

	precomputeStuff();
	precomputeLogKFact();

	getMatrixSize(argv[1]);

	allocateMemory();

	readFullOMatrix(argv[1], &O, true);
	readFullTMatrix(argv[2], &T);

	setBoundaries();
	if (argc == 3 || atoi(argv[3]) == 0)
		evaluateFullLikelihood();
	// else
	// 	evaluateFullLikelihoodSkipBoundary();
}
