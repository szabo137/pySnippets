#include <math.h>

double integF1(int n,double *x, void *userData){
    double c = *(double *)userData ;
    return cos(x[0])*(cos(x[0]/2.0/c)*cos(x[0]/2.0/c));
    }

double integF2(int n,double *x, void *userData){
    double c = *(double *)userData ;
    return sin(x[0])*(cos(x[0]/2.0/c)*cos(x[0]/2.0/c)) ;
    }

double integF3(int n,double *x, void *userData){
    double c = *(double *)userData ;
    return (cos(x[0])*(cos(x[0]/2.0/c)*cos(x[0]/2.0/c)))*(cos(x[0])*(cos(x[0]/2.0/c)*cos(x[0]/2.0/c)));
    }

double integF4(int n,double *x, void *userData){
    double c = *(double *)userData ;
    return (sin(x[0])*(cos(x[0]/2.0/c)*cos(x[0]/2.0/c)))*sin(x[0])*(cos(x[0]/2.0/c)*cos(x[0]/2.0/c)) ;
    }
