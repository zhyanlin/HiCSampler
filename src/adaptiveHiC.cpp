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

	if (option_gibbs)
	{ // gibbs is not implemented yet
		exit(1);
	}
	precomputeStuff();
	precomputeLogKFact();

	getMatrixSize(argv[1]);

	allocateMemory();

	readFullOMatrix(argv[1], &O, true);

	//read annotation file
	//file format: idx	chr	startbp	endbp	annoLen	annoGC annoMap
	if (strlen(option_annotationFile)>0){
		std::ifstream annotationFile(option_annotationFile);
		//0       chr20   31060000        31065000        3560    0.5940  1.0000
	int idx, startbp,endbp;
	float gc,len,map;
	string chr;
	while (annotationFile >> idx >> chr >> startbp >> endbp >> len >> gc >> map)
	{
	 anno.gc[idx]=gc;
	 anno.len[idx]=len/(endbp-startbp);
	 anno.map[idx]=map;
	}
	}


	// cout << "start setting boundary\n";
	// setBoundaries();
	// cout << "setting boundary done\n";

	if (method_fixed)
	{
		computeTMatrix_fixed();
	}
	if (method_kde)
	{
		KDE kde;
		kde.initialize();
		kde.runKDE();
		cout << "KDE done\n";
	}

	if (method_mrf)
	{
		cout << "MRF...\n";
		MRF mrf;
		mrf.initialize();
		mrf.runMRF();
		cout << "MRF done\n";
	}

	if (option_mcmc)
	{
		/* gamma potential based on kde
		if(strlen(option_prior_file)==0){
			//compute prior from T_kde 
			KDE kde;
			kde.initialize();
			kde.runKDE();
			double mean,s,k,theta;
			double eps=1e-10;
			int minSamples=50;//minimum number of values requried to compute prior
			int maxGenomicDistance=T.endRow()-T.startRow();
			double sumOfTs[maxGenomicDistance],sumOfLogTs[maxGenomicDistance];
			double sumOfT,sumOfLogT;
			int nSample;
			int nSamples[maxGenomicDistance];
			int offSet; int direction;
			for(int i=0;i<maxGenomicDistance;i++){
				sumOfTs[i]=0;
				sumOfLogTs[i]=0;
				nSamples[i]=0;
			}
			for (int i = T.startRow(); i < T.endRow(); i++)
			{
				for (int j = T.startCol(i); j < T.endCol(i); j++)
				{	
					if(T.get(i,j)>1e-4){
					sumOfTs[j-i]+=T.get(i,j);
					sumOfLogTs[j-i]+=log(T.get(i,j));
					nSamples[j-i]+=1;
					}		
				}
			}

			for (int i = 0; i < maxGenomicDistance; i++)
			{
				nSample=nSamples[i];
				sumOfT=sumOfTs[i];
				sumOfLogT=sumOfLogTs[i];
				offSet=1;
				direction=-1;
				while(nSample<50){
					if(i+offSet*direction>0 && i+offSet*direction<maxGenomicDistance){
						nSample+=nSamples[i+offSet*direction];
						sumOfT+=sumOfTs[i+offSet*direction];
						sumOfLogT+=sumOfLogTs[i+offSet*direction];
					}
					direction*=-1;
					if(direction==-1)
						offSet++;
				}
				mean = sumOfT/nSample;
				s=log(mean)-sumOfLogT/nSample;
		  		k=(3-s+sqrt(pow(s-3,2)+24*s))/(12*s+eps);
		  		theta=sumOfT/((k+eps)*nSample);

				piror_k[i] = k;
				piror_theta[i] = theta;
			}
		}*/
		cout << "MCMC...\n";

		adaptMCMC mcmc;
		mcmc.initialize(argv[2]);
		mcmc.burnIn();

		mcmc.mainIter();

		mcmc.outputSparseMatrixGZ(argv[2]);
		cout << "MCMC done\n";
	}
	else
		outputSparseMatrixGZ(argv[2]);
	if (option_testMatrix[0])
		evaluateAgainstTest();
}
