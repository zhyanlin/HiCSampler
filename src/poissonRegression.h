#ifndef POISSONREGRESSION_H_
#define POISSONREGRESSION_H_
#include <math.h>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
using namespace std;

struct GLMData
{
    int numOfFeatures;
    int offsets;
    int numOfInstances;
    double *y, *hicbias;
    double **X;
};
extern GLMData glmTraindata;
class poissonGLM
{
public:
    poissonGLM(const int features, const int numOfInstances, const int offsets = 0)
    {
        this->numOfFeatures = features;
        this->numOfOffsets = offsets;
        this->numOfInstances = numOfInstances;

        glmTraindata.numOfFeatures = features;
        glmTraindata.offsets = offsets;
        glmTraindata.numOfInstances = numOfInstances;

        this->w = new double[this->numOfFeatures + 1]; // w[-1] is the bias in linear regression
        this->wNew = new double[this->numOfFeatures + 1]; // w[-1] is the bias in linear regression
        this->y = new double[numOfInstances];
        this->X = new double *[numOfInstances];

        glmTraindata.y = new double[numOfInstances];
        glmTraindata.hicbias = new double[numOfInstances];
        glmTraindata.X = new double *[numOfInstances];
        for (int i = 0; i < numOfInstances; i++)
            glmTraindata.X[i] = new double[this->numOfFeatures + 1];

        for (int i = 0; i < this->numOfFeatures + 1; i++)
            this->w[i] = (float)rand() / RAND_MAX * 2 - 1;
        for (int i = 0; i < this->numOfOffsets; i++)
            this->w[this->numOfFeatures - 1 - i] = 1;

        this->mean = new double[this->numOfFeatures];
        this->std = new double[this->numOfFeatures];
        for (int i = 0; i < this->numOfFeatures; i++)
        {
            this->mean[i] = 0;
            this->std[i] = 0;
        }
    };
    int lbfgsfit(double **X, double *y,double *Xbias);
    int fit(double **X, double *y)
    {

        double lr = 1;
        double ll_new, ll_old;
        double wx;
        int inneriter;
        double *pd;
        pd = new double[this->numOfFeatures + 1];
        // cout<<numOfInstances<<"x"<<numOfFeatures<<endl;
        //standardization, get feature mean and std
        for (int i = 0; i < this->numOfInstances; i++)
            for (int j = 0; j < this->numOfFeatures; j++)
            {
                mean[j] += X[i][j];
                std[j] += pow(X[i][j], 2);
            }
        for (int i = 0; i < this->numOfFeatures; i++)
        {
            mean[i] /= this->numOfInstances;
            std[i] /= this->numOfInstances;
            std[i] = sqrt(std[i] - pow(mean[i], 2));
        }

        //apply the standardization
        for (int i = 0; i < this->numOfInstances; i++)
        {
            glmTraindata.y[i] = y[i];
            for (int j = 0; j < this->numOfFeatures; j++)
                glmTraindata.X[i][j] = (X[i][j] - mean[j]) / std[j];
            glmTraindata.X[i][this->numOfFeatures] = 1; // w[-1] is the bias in linear regression
        }
        //start fitting through gradient descent
        ll_old = 0;
        for (int i = 0; i < this->numOfInstances; i++)
        {
            wx = 0;
            for (int j = 0; j < this->numOfFeatures + 1; j++)
            {
                wx += w[j] * glmTraindata.X[i][j];
            }

            ll_old += glmTraindata.y[i] * wx - exp(wx);
        }

        //============
        ll_new = ll_old - 1;
        while (ll_new < ll_old)
        {
            // cout<<"before main iter lr="<<lr<<endl;
            lr /= 2.0;
            for (int i = 0; i < this->numOfFeatures + 1; i++)
                pd[i] = 0;

            for (int i = 0; i < this->numOfInstances; i++)
            {
                wx = 0;
                for (int j = 0; j < this->numOfFeatures + 1; j++)
                    wx += w[j] * glmTraindata.X[i][j];
                for (int j = 0; j < this->numOfFeatures + 1; j++)
                    pd[j] += glmTraindata.X[i][j] * exp(wx) - glmTraindata.y[i] * glmTraindata.X[i][j];
            }
            for (int i = 0; i < this->numOfFeatures + 1; i++)
            {
                if (i < this->numOfFeatures - this->numOfOffsets || i == this->numOfFeatures)
                {
                    wNew[i] = w[i] - lr * pd[i];
                    continue;
                }
                else
                    wNew[i] = 1;
            }

            ll_new = 0;
            for (int i = 0; i < this->numOfInstances; i++)
            {
                wx = 0;
                for (int j = 0; j < this->numOfFeatures + 1; j++)
                    wx += wNew[j] * glmTraindata.X[i][j];
                ll_new += glmTraindata.y[i] * wx - exp(wx);
            }
        } //end of adapting learning rate
          //============
        // lr=0.001;
        ll_new = ll_old - 1;
        while (fabs(ll_new - ll_old) > 1e-5)
        {   cout<<"ll_new "<<ll_new<<endl;

            ll_old = ll_new;
            //apply gradient descent
            for (int i = 0; i < this->numOfFeatures + 1; i++)
                pd[i] = 0;

            for (int i = 0; i < this->numOfInstances; i++)
            {
                wx = 0;
                for (int j = 0; j < this->numOfFeatures + 1; j++)
                    wx += w[j] * glmTraindata.X[i][j];
                for (int j = 0; j < this->numOfFeatures + 1; j++)
                    pd[j] += glmTraindata.X[i][j] * exp(wx) - glmTraindata.y[i] * glmTraindata.X[i][j];
            }
            for (int i = 0; i < this->numOfFeatures + 1; i++)
            {
                if (i < this->numOfFeatures - this->numOfOffsets || i == this->numOfFeatures)
                {
                    wNew[i] = w[i] - lr * pd[i];
                    // cout<<"wNew["<<i<<"] "<<wNew[i]<<" = "<<w[i]<<" - "<<lr<<"*"<<pd[i]<<endl;
                    continue;
                }
                else
                    wNew[i] = 1;
            }

            ll_new = 0;
            for (int i = 0; i < this->numOfInstances; i++)
            {
                wx = 0;
                for (int j = 0; j < this->numOfFeatures + 1; j++)
                    wx += wNew[j] * glmTraindata.X[i][j];
                ll_new += glmTraindata.y[i] * wx - exp(wx);
            }

            //end of gradient descent
            inneriter = 0;
            while (ll_new < ll_old && inneriter < 10)
            {
                cout << "lr=" << lr << endl;
                cout << ll_new << endl;
                inneriter++;
                lr /= 2.0;
                for (int i = 0; i < this->numOfFeatures + 1; i++)
                    pd[i] = 0;

                for (int i = 0; i < this->numOfInstances; i++)
                {
                    wx = 0;
                    for (int j = 0; j < this->numOfFeatures + 1; j++)
                        wx += w[j] * glmTraindata.X[i][j];
                    for (int j = 0; j < this->numOfFeatures + 1; j++)
                        pd[j] += glmTraindata.X[i][j] * exp(wx) - glmTraindata.y[i] * glmTraindata.X[i][j];
                }
                for (int i = 0; i < this->numOfFeatures + 1; i++)
                {
                    if (i < this->numOfFeatures - this->numOfOffsets || i == this->numOfFeatures)
                    {
                        wNew[i] = w[i] - lr * pd[i];
                        continue;
                    }
                    else
                        wNew[i] = 1;
                }

                ll_new = 0;
                for (int i = 0; i < this->numOfInstances; i++)
                {
                    wx = 0;
                    for (int j = 0; j < this->numOfFeatures + 1; j++)
                        wx += wNew[j] * glmTraindata.X[i][j];
                    ll_new += glmTraindata.y[i] * wx - exp(wx);
                }
            } //end of adapting learning rate
            cout << ll_old << " -> " << ll_new << " , " << ll_new - ll_old << "(lr=" << lr << ")" << endl;
            for (int i = 0; i < this->numOfFeatures + 1; i++)
                w[i] = wNew[i];
        }
        cout << "ll_new " << ll_new << endl;
        for (int i = 0; i < this->numOfFeatures + 1; i++)
        cout<<"w["<<i<<"]="<<w[i]<<endl;
    };
    double predict(double *x)
    {
        double val = 0;
        for (int i = 0; i < this->numOfFeatures; i++){
            // cout<<"w[i] * ((x["<<i<<"] - this->mean["<<i<<"]) / this->std["<<i<<"])"<<endl;
            // cout<<" > "<<w[i]<<" *"<< "(("<<x[i]<<" -"<< this->mean[i]<<") /"<< this->std[i]<<")"<<endl;
            val += w[i] * ((x[i] - this->mean[i]) / this->std[i]);
        
        }
        val += w[this->numOfFeatures]; //GLM bias
        // cout<<"w[this->numOfFeatures] "<<w[this->numOfFeatures]<<endl;
        return exp(val);
    };

private:
    int numOfFeatures, numOfOffsets, numOfInstances;
    double *w, *mean, *std, *aggInstance, *wNew;
    double **X, *y;
};

#endif