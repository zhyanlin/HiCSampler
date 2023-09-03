#include "readOptions.h"
#include "utils.h"
#include "KDE.h"
#include "MRF.h"
#include "adaptMCMC.h"
#include "contactMap.h"
#include "optimize.h"
#include <fstream>
using namespace std;

int main(int argc, char *argv[])
{
	srand48(time(NULL));
	if (argc < 4)
	{
		fprintf(stderr,
				"adaptiveHiC <matrix_file_with_row_count> outputfile <options>\n");
		exit(1);
	}

	readOptions(argc, argv);

	FILE *temp = fopen(argv[2], "w");
	if (temp == NULL)
	{
		fprintf(stderr, "Error: Can't open output file: %s\n", argv[2]);
		exit(1);
	}
	fclose(temp);

	precomputeStuff();
	precomputeLogKFact();

	getMatrixSize(argv[1]);

	allocateMemory();

	readFullOMatrix(argv[1], &O, true);

	setBoundaries();

	computeTMatrix_fixed();

	cout << "MCMC...\n";

	adaptMCMC mcmc;
	mcmc.initialize(argv[2]);
	mcmc.burnIn();

	mcmc.mainIter();

	mcmc.outputSparseMatrixGZ(argv[2]);
	cout << "MCMC done\n";
}
