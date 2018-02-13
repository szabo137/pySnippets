#include <math.h>

double g(double x){
    return exp(x)*sin(3*x);
    }

double f(int n, double *x, void *userData){
    double c= *(double *)userData;
    return c + g(x[0]) - x[1]*x[2];
    }
