#include "test-functions.h"
#include "./../core/constants.h"
#include <math.h>
#include <iostream>

double c2_discontinuity(double x, double y, double param, double* domain){
    if(x <= param)
        return (x+domain[0])*(x+domain[0]);
    else if(x > param)
        return (2*param*(x+domain[0]) - param*param);
}

double c1_discontinuity(double x, double y, double param, double* domain){
    if( x <= param)
        return (x)*(x);
    else
        return ((param*param) / (1 - param)) * (1 - (x));
}

double constant_linear_constant(double x, double y, double* domain){

//    x = x + y*1.373;
    if(fabs(x) <= 0.3)
        return  0.7 - x;
    else if(x < -0.3)
        return 1.0;
    else if(x > 0.3)
        return 0.4;
}

double rotated_constant_linear_constant(double x, double y, double* domain){

    x = x + y*1.373;
    if(fabs(x) <= 0.3)
        return  0.7 - x;
    else if(x < -0.3)
        return 1.0;
    else if(x > 0.3)
        return 0.4;
}

double constant_linear(double x, double y, double* domain){
    if(x < -0.29)
        return 1.0;
    else
        return 0.71 - x;
}

double linear_linear(double x, double y, double* domain){
    if(x < -0.313)
        return 1.313 + x;
    else
        return 0.687 - x;
}

double triangular_distribution(double x, double xPeakPos, double* domain){

    double a = domain[0];
    double b = domain[2];
    double c = xPeakPos;

    if( x < a) return 0;
    if(x >= a && x < c)
        return (2*(x - a)) / ((b-a)*(c-a));
    if(x == c)
        return 2 / (b-a);
    if(x > c && x <= b)
        return (2*(b-x))/((b-a)*(b-c));
    if(x > b)
        return 0;
}

double linear(double x, double y, double* domain){
    return 1.414 - 1.373*x;
}

double plane(double x, double y, double* domain){
    return 1.414 + x + 1.373*y;
}

double quadratic(double x, double y, double* domain){
    return 1.414 + x*x + 1.373*y*y;
}

