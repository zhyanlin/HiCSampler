#include "readOptions.h"
int n_threads = 4;
int option_maxFragDistance = 0;
int k = 100;
int option_mcmcMaxIter = 100000;
char option_method[100];
bool method_mrf = false;
int Tr = 0;
int Tc = 0;
bool method_fixed = false;
bool method_kde = false;
bool option_lognormal = true;
bool option_median = true;
double option_kdeBandwidth = -1;
double option_kdeMaxBandwidth = 25;
double option_minSigma = 0.001;
double option_sigmaMultiplier = 0.04;
bool option_kdeInitialization = false;
int option_diagBand = 1000;
int option_kdeMinCount = 50;
double option_boundaryDensity = 100.0;
double option_boundaryDenisty = 1.0;
double option_sigma = 10;
double option_averageFragmentSize = 4096;
bool option_inputFullMatrix = false;
bool option_inputSparseMatrix = true;
bool option_noBias = false;
double option_pseudocount = 0.01;
bool option_mcmc = false;
bool option_gibbs = false;
double option_boundaryKS = 2.0;
bool option_diagRenormalization = true;
bool option_outputNormalized = false;
double option_minOutput = 0.00001;
int option_mrfMaxIter = 10;
int option_firstRow = -1;
int option_firstCol = -1;
int option_lastRow = -1;
int option_lastCol = -1;
char option_initializationMatrix[1000];
char option_testMatrix[1000];

char option_annotationFile[1000];

char option_bias[1000];
char option_prior_file[1000];
char option_boundaryOutput[1000];
bool option_outputSymmetric = false;
int bandSize = 9999999;
int stepSize = 500;
bool outputSamples=true;
bool option_diffusion = false;
#include <iostream>

void readOptions(int argc, char *argv[])
{
  for (int i = 3; i < argc; i++)
  {

    if (strstr(argv[i], "-averageFragmentSize="))
    {
      int x = sscanf(argv[i], "-averageFragmentSize=%lf",
                     &option_averageFragmentSize);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }

    if (strstr(argv[i], "-method="))
    {
      int x = sscanf(argv[i], "-method=%s", option_method);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      if (!strcmp(option_method, "mrf"))
      {
        method_mrf = true;
        continue;
      }
      if (!strcmp(option_method, "kde"))
      {
        method_kde = true;
        method_mrf = false;
        continue;
      }
      if (!strcmp(option_method, "fixed"))
      {
        method_fixed = true;
        method_mrf = false;
        continue;
      }
      fprintf(stderr,
              "Unknown method: %s\nChoices are: mrf, kde, fixed\n",
              option_method);
    }
    if (strstr(argv[i], "-initializationMatrix"))
    {
      sscanf(argv[i], "-initializationMatrix=%s",
             option_initializationMatrix);
      continue;
    }

    if (strstr(argv[i], "-kdeBandwidth="))
    {
      int x = sscanf(argv[i], "-kdeBandwidth=%lf", &option_kdeBandwidth);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }
    if (strstr(argv[i], "-kdeMaxBandwidth="))
    {
      int x = sscanf(argv[i], "-kdeMaxBandwidth=%lf",
                     &option_kdeMaxBandwidth);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }
    if (strstr(argv[i], "-kdeMinCount="))
    {
      int x = sscanf(argv[i], "-kdeMinCount=%d", &option_kdeMinCount);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }
    if (strstr(argv[i], "-boundaryDensity"))
    {
      int x = sscanf(argv[i], "-boundaryDensity=%lf",
                     &option_boundaryDensity);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }
    if (strstr(argv[i], "-sigma="))
    {
      int x = sscanf(argv[i], "-sigma=%lf", &option_sigma);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }
    if (strstr(argv[i], "-boundaryKS="))
    {
      int x = sscanf(argv[i], "-boundaryKS=%lf", &option_boundaryKS);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }

    if (strstr(argv[i], "-k="))
    {
      int x = sscanf(argv[i], "-k=%d", &k);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }

    if (strstr(argv[i], "-Tr="))
    {
      int x = sscanf(argv[i], "-Tr=%d", &Tr);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }

    if (strstr(argv[i], "-Tc="))
    {
      int x = sscanf(argv[i], "-Tc=%d", &Tc);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      }
      continue;
    }

    if (!strcmp(argv[i], "-mcmc"))
    {
      option_mcmc = true;
      continue;
    }
    if (!strcmp(argv[i], "-new"))
    {
      option_diffusion = true;
      continue;
    }
    if (!strcmp(argv[i], "-gibbs"))
    {
      option_gibbs = true;
      continue;
    }
    if (!strcmp(argv[i], "-kdeInitialization"))
    {
      option_kdeInitialization = true;
      continue;
    }

    if (!strcmp(argv[i], "-inputFullMatrix"))
    {
      option_inputFullMatrix = true;
      option_inputSparseMatrix = false;
      continue;
    }
    if (!strcmp(argv[i], "-inputSparseMatrix"))
    {
      option_inputSparseMatrix = true;
      option_inputFullMatrix = false;
      continue;
    }

    if (!strcmp(argv[i], "-outputNormalized"))
    {
      option_outputNormalized = true;
      continue;
    }
    if (!strcmp(argv[i], "-noBias"))
    {
      option_noBias = true;
      continue;
    }
    if (strstr(argv[i], "-minOutput="))
    {
      int x = sscanf(argv[i], "-minOutput=%lf", &option_minOutput);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }
    if (strstr(argv[i], "-pseudocount="))
    {
      int x = sscanf(argv[i], "-pseudocount=%lf", &option_pseudocount);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }

    if (strstr(argv[i], "-mcmcIter="))
    {
      int x = sscanf(argv[i], "-mcmcIter=%d", &option_mcmcMaxIter);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }

    if (strstr(argv[i], "-mrfMaxIter="))
    {
      int x = sscanf(argv[i], "-mrfMaxIter=%d", &option_mrfMaxIter);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }

    if (strstr(argv[i], "-threads="))
    {
      int x = sscanf(argv[i], "-threads=%d", &n_threads);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      std::cout << "n_threads=" << n_threads << std::endl;
      continue;
    }

    if (strstr(argv[i], "-bandSize="))
    {
      int x = sscanf(argv[i], "-bandSize=%d", &bandSize);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      std::cout << "bandSize=" << bandSize << std::endl;
      continue;
    }

    if (strstr(argv[i], "-stepSize="))
    {
      int x = sscanf(argv[i], "-stepSize=%d", &stepSize);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      std::cout << "stepSize=" << stepSize << std::endl;
      continue;
    }

    if (strstr(argv[i], "-firstRow="))
    {
      int x = sscanf(argv[i], "-firstRow=%d", &option_firstRow);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }
    if (strstr(argv[i], "-lastRow="))
    {
      int x = sscanf(argv[i], "-lastRow=%d", &option_lastRow);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }

    if (strstr(argv[i], "-firstColumn="))
    {
      int x = sscanf(argv[i], "-firstColumn=%d", &option_firstCol);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }
    if (strstr(argv[i], "-lastColumn="))
    {
      int x = sscanf(argv[i], "-lastColumn=%d", &option_lastCol);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }
    if (strstr(argv[i], "-testMatrix="))
    {
      int x = sscanf(argv[i], "-testMatrix=%s", option_testMatrix);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }

    if (strstr(argv[i], "-a="))
    {
      int x = sscanf(argv[i], "-a=%s", option_annotationFile);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }


    if (strstr(argv[i], "-bias="))
    {
      int x = sscanf(argv[i], "-bias=%s", option_bias);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }


        if (strstr(argv[i], "-prior="))
    {
      int x = sscanf(argv[i], "-prior=%s", option_prior_file);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }

    if (strstr(argv[i], "-boundaryOutput="))
    {
      int x = sscanf(argv[i], "-boundaryOutput=%s",
                     option_boundaryOutput);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }
    if (strstr(argv[i], "-outputSymmetric"))
    {
      option_outputSymmetric = true;
      continue;
    }
    if (strstr(argv[i], "-lognormal"))
    {
      option_lognormal = true;
      continue;
    }
    if (strstr(argv[i], "-median"))
    {
      option_median = true;
      continue;
    }
    if (strstr(argv[i], "-minSigma="))
    {
      int x = sscanf(argv[i], "-minSigma=%lf", &option_minSigma);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }
    if (strstr(argv[i], "-sigmaMultiplier="))
    {
      int x = sscanf(argv[i], "-sigmaMultiplier=%lf",
                     &option_sigmaMultiplier);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }
    if (strstr(argv[i], "-diagBand="))
    {
      int x = sscanf(argv[i], "-diagBand=%d", &option_diagBand);
      if (x != 1)
      {
        fprintf(stderr, "Error: Can't read %s\n", argv[i]);
      };
      continue;
    }

    fprintf(stderr, "Unknown option: %s\n", argv[i]);
    exit(1);
  }

  printf("option -method= %s\n", option_method);
  printf("option -kdeMinCount = %d\n", option_kdeMinCount);

  printf("option -kdeBandwidth = %lf\n", option_kdeBandwidth);
  printf("option -kdeMaxBandwidth = %lf\n", option_kdeMaxBandwidth);
  printf("option -boundaryDensity = %lf\n", option_boundaryDensity);
  printf("option -kdeInitialization = %d\n", (int)option_kdeInitialization);
  printf("option -sigma = %lf\n", option_sigma);
  printf("option -averageFragmentSize = %lf\n", option_averageFragmentSize);
  printf("option -inputFullMatrix = %d\n", (int)option_inputFullMatrix);
  printf("option -inputSparseMatrix= %d\n", (int)option_inputSparseMatrix);
  printf("option -minOutput = %lf\n", option_minOutput);
  printf("option -pseudocount = %lf\n", option_pseudocount);
  printf("option -mrfMaxIter = %d\n", option_mrfMaxIter);
  printf("option -firstRow = %d\n", option_firstRow);
  printf("option -lastRow = %d\n", option_lastRow);
  printf("option -firstCol = %d\n", option_firstCol);
  printf("option -lastCol = %d\n", option_lastCol);
  printf("option -initializationMatrix = %s\n", option_initializationMatrix);
  printf("option -minSigma = %f\n", option_minSigma);
  printf("option -sigmaMultiplier = %f\n", option_sigmaMultiplier);

  printf("option -testMatrix = %s\n", option_testMatrix);
  printf("option -gibbs = %d\n", (int)option_gibbs);
  printf("option -mcmc = %d\n", (int)option_mcmc);

  printf("option -lognormal = %d\n", (int)option_lognormal);
  printf("option -median = %d\n", (int)option_median);
  printf("option -minSigma = %f\n", option_minSigma);
  printf("option -boundaryKS = %f\n", option_boundaryKS);
};
