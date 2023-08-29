#include "../readOptions.h"
#include "../utils.h"
#include "../contactMap.h"
#include <stdlib.h>
#include <gsl/gsl_randist.h>
using namespace std;

double evaluateLikelihood(const float& currentT,const  int& obs,const  float& medNei,const  float& bias) {
  float transformedMedNeighbor;

  transformedMedNeighbor=log(medNei+1);
  // no +1 here?
  float transformedCurrentTbias;
  transformedCurrentTbias=log((currentT+0.00001)*bias);

  float localSigma=max(option_minSigma, sqrt(medNei*option_sigmaMultiplier));
  float NORM_RES_over_sigma=NORM_RES/localSigma;

  float correctedp=currentT*bias;

  // log of Poisson distribution
  float loglike1= -correctedp + obs  * (transformedCurrentTbias) - logKFact[obs];

  float transformedp;
  transformedp=log(currentT+1);

  float loglike2=0;
  int iii=(transformedp-transformedMedNeighbor)*NORM_RES_over_sigma + fullRange;
  cout<<"("<<transformedp<<" - "<<transformedMedNeighbor<<")*"<<NORM_RES_over_sigma<<"+"<<fullRange<<" = "<<(transformedp-transformedMedNeighbor)*NORM_RES_over_sigma + fullRange;
  cout<<endl;
  cout<<"< "<<floatFullRange<<"?\n";
  if (iii<0) iii=0;
  if (iii>=floatFullRange) iii=floatFullRange-1;
  loglike2= precomputedNormLogPDF[iii];
  cout<<iii<<endl;
  cout<<loglike1<<" "<<loglike2<<" "<<loglike1+loglike2<<endl;
  return loglike1+loglike2;
}

int main(int argc, char *argv[])
{

	srand48(time(NULL));
	if (argc != 5)
	{
		fprintf(stderr,
				"decomp_ll Obs IF median Bias\n");
		exit(1);
	}


	precomputeStuff();
	precomputeLogKFact();

	double medNei,bias,IF;
	int obs;
	medNei=atof(argv[3]);
	IF=atof(argv[2]);
	bias=atof(argv[4]);
	obs=atoi(argv[1]);
	evaluateLikelihood(IF,obs,medNei,bias);


	return 0;


}
