#include "curvefit.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
//y=a*e^(bx)+c

struct FITEXPDATA
{
    double *x;
    double *y;
    int n;
} fitexpData;

void exp_df(const gsl_vector *w, void *params,
            gsl_vector *df)
{
    double dlda, dldb, dldc;
    dlda = 0;
    dldb = 0;
    dldc = 0;
    double a, b, c;
    a = gsl_vector_get(w, 0);
    b = gsl_vector_get(w, 1);
    c = gsl_vector_get(w, 2);
    for (int i = 0; i < fitexpData.n; i++)
    {
        dlda += 2 * a * exp(2 * b * fitexpData.x[i]) + c * c + fitexpData.y[i] * fitexpData.y[i] + 2 * (c - fitexpData.y[i]) * exp(b * fitexpData.x[i]) - 2 * fitexpData.y[i] * c;

        dldb += 2 * fitexpData.x[i] * a * a * exp(2 * b * fitexpData.x[i]) + c * c + fitexpData.y[i] * fitexpData.y[i] +
                2 * a * fitexpData.x[i] * (c - fitexpData.y[i]) * exp(b * fitexpData.x[i]) - 2 * c * fitexpData.y[i];

        // cout << "dldc[" << i << "]= " << dldc << endl;
        dldc += a * a * exp(2 * b * fitexpData.x[i]) + 2 * c + fitexpData.y[i] * fitexpData.y[i] + 2 * a * exp(b * fitexpData.x[i]) -
                2 * a * fitexpData.y[i] * exp(b * fitexpData.x[i]) - 2 * fitexpData.y[i];
    }
    // cout << "dldc= " << dldc << endl;
    // cout << a << " " << b << " " << c << endl;
    gsl_vector_set(df, 0, dlda);
    gsl_vector_set(df, 1, dldb);
    gsl_vector_set(df, 2, dldc);
};

double exp_f(const gsl_vector *w, void *params)
{
    double f = 0;

    for (int i = 0; i < fitexpData.n; i++)
    {
        f += pow(gsl_vector_get(w, 0) * exp(gsl_vector_get(w, 1) * fitexpData.x[i]) + gsl_vector_get(w, 2) - fitexpData.y[i], 2);
    }
    cout << "(exp) f=" << f << endl;
    return f;
};

void exp_fdf(const gsl_vector *x, void *params,
             double *f, gsl_vector *df)
{

    *f = exp_f(x, params);
    exp_df(x, params, df);
};

//y=a*e^(bx)+c
void fitExp_bfgs(int n, const double *x, const double *Y, double &a, double &b, double &c)
{

    fitexpData.x = new double[n];
    fitexpData.y = new double[n];
    fitexpData.n = n;
    for (int i = 0; i < n; i++)
    {
        fitexpData.x[i] = x[i];
        fitexpData.y[i] = Y[i];
    }
    double par[1] = {0};

    gsl_vector *w;
    gsl_multimin_function_fdf my_func;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    // function
    my_func.n = 3;
    my_func.f = exp_f;
    my_func.df = exp_df;
    my_func.fdf = exp_fdf;
    my_func.params = par;
    /* Starting point, x = (5,7) */
    w = gsl_vector_alloc(3);

    gsl_vector_set(w, 0, a);
    gsl_vector_set(w, 1, b);
    gsl_vector_set(w, 2, c);

    T = gsl_multimin_fdfminimizer_conjugate_pr;
    s = gsl_multimin_fdfminimizer_alloc(T, 3);

    gsl_multimin_fdfminimizer_set(s, &my_func, w, 1, 1e-1);
    int iter = 0;
    int status;

    do
    {
        cout << "iter " << iter << endl;
        cout << "s->x " << gsl_vector_get(s->x, 0) << endl;
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

               if (status)
                 break;

        status = gsl_multimin_test_gradient(s->gradient, 1e-10);

    } while (status == GSL_CONTINUE && iter < 1000);

    a = gsl_vector_get(s->x, 0);
    b = gsl_vector_get(s->x, 1);
    c = gsl_vector_get(s->x, 2);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(w);
    delete[] fitexpData.x;
    delete[] fitexpData.y;
};
double predExp(double x, double a, double b, double c)
{
    return a * exp(b * x) + c;
}
void fitExpFinetune(int n, const double *x, const double *y, double &a, double &b, double &c, bool tune_a, bool tune_b, bool tune_c)
{
    double lra, lrb, lrc;
    lra = 0;
    lrb = 0;
    lrc = 0;
    if (tune_a)
        lra = 1e-2;
    if (tune_b)
        lrb = 1e-2;
    if (tune_c)
        lrc = 1e-2;

    double eps = 1e-12;
    double dlda, dldb, dldc;
    double eval_old, eval_new, eval_tmp;
    double a_new, b_new, c_new;
    eval_old = 0;
    eval_new = 0;
    for (int i = 0; i < n; i++)
    {
        eval_old += pow(a * exp(b * x[i]) + c - y[i], 2);
    }

    dlda = 0;
    dldb = 0;
    dldc = 0;
    for (int i = 0; i < n; i++)
    {
        dlda += 2 * a * exp(2 * b * x[i]) + c * c + y[i] * y[i] + 2 * (c - y[i]) * exp(b * x[i]) - 2 * y[i] * c;

        dldb += 2 * x[i] * a * a * exp(2 * b * x[i]) + c * c + y[i] * y[i] + 2 * a * x[i] * (c - y[i]) * exp(b * x[i]) - 2 * c * y[i];

        dldc += a * a * exp(2 * b * x[i]) + 2 * c + y[i] * y[i] + 2 * a * exp(b * x[i]) - 2 * a * y[i] * exp(b * x[i]) - 2 * y[i];
    }
    a_new = a - lra * dlda;
    b_new = b - lrb * dldb;
    c_new = c - lrc * dldc;

    for (int i = 0; i < n; i++)
    {
        eval_new += pow(a_new * exp(b_new * x[i]) + c_new - y[i], 2);
    }

    while (eval_new > eval_old)
    {
        lra /= 2;
        lrb /= 2;
        lrc /= 2;
        eval_new = 0;
        dlda = 0;
        dldb = 0;
        dldc = 0;
        for (int i = 0; i < n; i++)
        {
            dlda += 2 * a * exp(2 * b * x[i]) + c * c + y[i] * y[i] + 2 * (c - y[i]) * exp(b * x[i]) - 2 * y[i] * c;

            dldb += 2 * x[i] * a * a * exp(2 * b * x[i]) + c * c + y[i] * y[i] + 2 * a * x[i] * (c - y[i]) * exp(b * x[i]) - 2 * c * y[i];

            dldc += a * a * exp(2 * b * x[i]) + 2 * c + y[i] * y[i] + 2 * a * exp(b * x[i]) - 2 * a * y[i] * exp(b * x[i]) - 2 * y[i];
        }
        a_new = a - lra * dlda;
        b_new = b - lrb * dldb;
        c_new = c - lrc * dldc;

        for (int i = 0; i < n; i++)
        {
            eval_new += pow(a_new * exp(b_new * x[i]) + c_new - y[i], 2);
        }
    }
    // . optim lr search
    eval_new = 0;
    for (int i = 0; i < n; i++)
    {
        eval_new += pow(a * exp(b * x[i]) + c - y[i], 2);
    }

    // a_new=100;b_new=100;c_new=100;
    eval_new = eval_old;
    eval_old = 1e32;
    int iter = 0;
    while (fabs(eval_new - eval_old) > eps)
    {
        cout << "eval_new=" << eval_new << "; eval_old=" << eval_old << endl;
        iter++;
        eval_tmp = eval_old;
        eval_old = eval_new;
        dlda = 0;
        dldb = 0;
        dldc = 0;
        for (int i = 0; i < n; i++)
        {
            dlda += 2 * a * exp(2 * b * x[i]) + c * c + y[i] * y[i] + 2 * (c - y[i]) * exp(b * x[i]) - 2 * y[i] * c;

            dldb += 2 * x[i] * a * a * exp(2 * b * x[i]) + c * c + y[i] * y[i] + 2 * a * x[i] * (c - y[i]) * exp(b * x[i]) - 2 * c * y[i];

            dldc += a * a * exp(2 * b * x[i]) + 2 * c + y[i] * y[i] + 2 * a * exp(b * x[i]) - 2 * a * y[i] * exp(b * x[i]) - 2 * y[i];
        }

        a_new = a - lra * dlda;
        b_new = b - lrb * dldb;
        c_new = c - lrc * dldc;
        eval_new = 0;

        for (int i = 0; i < n; i++)
        {
            eval_new += pow(a_new * exp(b_new * x[i]) + c_new - y[i], 2);
        }

        a = a_new;
        b = b_new;
        c = c_new;

        if (eval_new > eval_old)
            break;
        else
        {
            a = a_new;
            b = b_new;
            c = c_new;
        }
    }
}

//y=a*e^(bx)+c
void fitExp(int n, const double *x, const double *Y, double &a, double &b, double &c)
{

    double *y;
    y = new double[n];

    double sxy, sy, sxxy, sylny, sxylny;
    double c_old, c_delta;
    double sse_old, sse_new;
    sxy = 0;
    sy = 0;
    sxxy = 0;
    sylny = 0;
    sxylny = 0;

    for (int i = 0; i < n; i++)
        y[i] = Y[i] - c;

    for (int i = 0; i < n; i++)
    {
        sxy += x[i] * y[i];
        sy += y[i];
        sxxy += x[i] * x[i] * y[i];
        sylny += y[i] * log(y[i]);
        sxylny += x[i] * y[i] * log(y[i]);
    }

    a = (sxxy * sylny - sxy * sxylny) / (sy * sxxy - sxy * sxy);
    b = (sy * sxylny - sxy * sylny) / (sy * sxxy - sxy * sxy);
    a = exp(a);

    sse_old = 0;
    for (int i = 0; i < n; i++)
        sse_old += pow(Y[i] - c - a * exp(b * x[i]), 2);

    c_old = c;

    int maxiter = 1000;

    for (int iter = 0; iter < maxiter; iter++)
    {
        cout << "iter " << iter << endl;

        sxy = 0;
        sy = 0;
        sxxy = 0;
        sylny = 0;
        sxylny = 0;

        for (int i = 0; i < n; i++)
            y[i] = fmax(1e-10, Y[i] - c);
        for (int i = 0; i < n; i++)
        {
            sxy += x[i] * y[i];
            sy += y[i];
            sxxy += x[i] * x[i] * y[i];
            sylny += y[i] * log(y[i]);
            sxylny += x[i] * y[i] * log(y[i]);
        }
        a = (sxxy * sylny - sxy * sxylny) / (sy * sxxy - sxy * sxy);
        b = (sy * sxylny - sxy * sylny) / (sy * sxxy - sxy * sxy);
        a = exp(a);

        sse_new = 0;
        for (int i = 0; i < n; i++)
            sse_new += pow(Y[i] - c - a * exp(b * x[i]), 2);

        if (sse_new > sse_old)
        {
            c = c_old;
            break;
        }
        else
        {
            sse_old = sse_new;

            c_old = c;
            c_delta = 0;
            for (int i = 0; i < n; i++)
                c_delta += y[i] - a * exp(b * x[i]);
            c -= c_delta / n;
        }
    }
    delete[] y;

    fitExp_bfgs(n, x, Y, a, b, c);
    // fitExpFinetune(n, x, Y, a, b, c, 1, 1, 1);
}
