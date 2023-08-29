#include "poissonRegression.h"


GLMData glmTraindata;

/* The gradient of f, df = (df/dx, df/dy). */
void mle_df(const gsl_vector *w, void *params,
            gsl_vector *df)
{
    double *pd;
    pd = new double[glmTraindata.numOfFeatures + 1];
    for (int i = 0; i < glmTraindata.numOfFeatures + 1; i++)
        pd[i] = 0;

    double wx;
    for (int i = 0; i < glmTraindata.numOfInstances; i++)
    {
        wx = 0;
        for (int j = 0; j < glmTraindata.numOfFeatures + 1; j++)
            wx += gsl_vector_get(w, j) * glmTraindata.X[i][j];
    
        for (int j = 0; j < glmTraindata.numOfFeatures + 1; j++)
            pd[j] += glmTraindata.X[i][j] * exp(wx) - glmTraindata.y[i] * glmTraindata.X[i][j];
    }
    //reset gradient to 0 for offset features
    for (int i = 0; i < glmTraindata.offsets; i++)
        pd[glmTraindata.numOfFeatures - 1 - i] = 0;

    for (int i = 0; i < glmTraindata.numOfFeatures + 1; i++)
        gsl_vector_set(df, i, pd[i]);

    delete pd;
};

double mle_f(const gsl_vector *w, void *params)
{
    // for (int j = 0; j < glmTraindata.numOfFeatures + 1; j++)
    //         cout<<"w["<<j<<"] "<<gsl_vector_get(w, j)<<endl; 
    double wx, f;
    f = 0;
    for (int i = 0; i < glmTraindata.numOfInstances; i++)
    {
        wx = 0;
        for (int j = 0; j < glmTraindata.numOfFeatures + 1; j++)
            wx += gsl_vector_get(w, j) * glmTraindata.X[i][j];
        f += glmTraindata.y[i] * wx - exp(wx);
    }
    // cout<<"ll="<<f<<endl;
    f=-f;
    
    return f;
};

/* Compute both f and df together. */
void mle_fdf(const gsl_vector *x, void *params,
             double *f, gsl_vector *df)
{

    *f = mle_f(x, params);
    mle_df(x, params, df);
    
};


int poissonGLM::lbfgsfit(double **X, double *y)
{
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
        glmTraindata.y[i]= y[i];
        for (int j = 0; j < this->numOfFeatures; j++)
            glmTraindata.X[i][j] = (X[i][j] - mean[j]) / std[j];
        glmTraindata.X[i][this->numOfFeatures] = 1;
    }

    size_t iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

  
    double par[1] = {0};

    gsl_vector *x;
    gsl_multimin_function_fdf my_func;

    // function
    my_func.n = this->numOfFeatures + 1;
    my_func.f = mle_f;
    my_func.df = mle_df;
    my_func.fdf = mle_fdf;
    my_func.params = par;

    /* Starting point, x = (5,7) */
    x = gsl_vector_alloc(this->numOfFeatures + 1);
    for (int i = 0; i < this->numOfFeatures + 1; i++)
        gsl_vector_set(x, i,this->w[i]);

    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T, this->numOfFeatures + 1);

    gsl_multimin_fdfminimizer_set(s, &my_func, x, 1, 1e-2);
   
    do
    {   
        cout<<"iter "<<iter<<endl;
        iter++;
        status=gsl_multimin_fdfminimizer_iterate(s);

        if (status)
            break;

        status = gsl_multimin_test_gradient(s->gradient, 1e-10);
 


    } while (status == GSL_CONTINUE &&  iter < 10000);

    

    for (int i = 0; i < this->numOfFeatures + 1; i++){
        this->w[i]=gsl_vector_get(s->x, i);
        cout<<"w["<<i<<"]="<<this->w[i]<<endl;
        }
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

    return 0;
};