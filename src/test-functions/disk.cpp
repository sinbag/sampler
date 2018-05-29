#include "test-functions.h"

#include <math.h>
#include <iostream>

double disk_points2d(double x, double y, double diskRadius, double centerx, double centery,
                     double minvalue, double maxvalue){
    double dx = x - centerx;
    double dy = y - centery;
    double radius = sqrt(dx*dx + dy*dy);
    if(radius < diskRadius)
        return maxvalue;
    else
        return minvalue;
}

///C1 discontinous function
double disk_gaussian_points2d(double x, double y, double x1, double y1, double xcenter,
                              double ycenter, double parameter){

    double gaussVal = gaussian_points2d(x1, y1, xcenter, ycenter, parameter);
    double gauss  = gaussian_points2d(x, y, xcenter, ycenter, parameter);
    if(gauss > gaussVal)
        return gaussVal;
    else{
        return gauss;
    }
}

double disk_gaussian_points1d(double x, double x1, double xcenter, double parameter){

    double gaussVal = gaussian_points1d(x1, xcenter, parameter);
    double gauss  = gaussian_points1d(x, xcenter, parameter);
    if(gauss > gaussVal)
        return gaussVal;
    else{
        return gauss;
    }
}
