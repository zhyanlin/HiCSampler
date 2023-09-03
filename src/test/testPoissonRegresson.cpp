#include "../poissonRegression.h"
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    //     id,num_awards,prog,math
    // 45,0,3,41

    int features,instances;
    instances=0;
    features=2;
    FILE *f = fopen("./poissonData.csv", "r");
    char line[1000];
    double *y, **X;
    y = new double[1000];
    X = new double *[1000];
    for (int i = 0; i < 1000; i++)
        X[i] = new double[10];
    fscanf(f, "%[^\n]\n", line); //skip header

    int id;
    double num_awards,prog,math;
    while (fscanf(f, "%[^\n]\n", line) != EOF)
    {
        sscanf(line, "%d,%lf,%lf,%lf", &id,&num_awards, &prog, &math);        
        y[instances]=num_awards;
        X[instances][0]=prog;
        X[instances][1]=math;
        instances++;
    }
    poissonGLM pglm(features,instances,0);
    cout<<"start fitting\n";
    pglm.fit(X,y);


    poissonGLM pglm2(features,instances,1);
    cout<<"start fitting with offset\n";
    pglm2.fit(X,y);
};