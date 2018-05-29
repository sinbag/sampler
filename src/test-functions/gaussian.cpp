#include "test-functions.h"

#include <math.h>
#include <iostream>

double gaussian_points2d(double x, double y, double cx, double cy, double parameter){
    double normalizationFactor = 1.0;//1.0/(2*PI*parameter*parameter);
    double dx =  x - cx;
    double dy = y - cy;
    double gauss = normalizationFactor *exp(-(dx*dx+dy*dy)/(2*parameter*parameter));
    return gauss;
}

double gaussian_points1d(double x, double xcenter, double parameter){
    double normalizationFactor = 1.0;
    double dx = x - xcenter;
    double gauss = normalizationFactor *exp(-(dx*dx)/(2*parameter*parameter));
    return gauss;
}

double gaussian_mixture(double x, double y, double cx, double cy, int nGaussians, double* bbox, std::string model){


    double maxRange = bbox[2] - bbox[0];
    double minX = bbox[0];
    double minY = bbox[1];

    double dx =  x - cx;
    double dy = y - cy;
    double mixture = 0.0;
    for(int k = 0; k < nGaussians; k++){

        if(model=="random"){
            dx = x - (drand48()*maxRange + minX);
            dy = y - (drand48()*maxRange + minY);
        }

        double parameter = drand48();
        double normalizationFactor = 1.0;//1.0/(2*PI*parameter*parameter);
        mixture += normalizationFactor *exp(-(dx*dx+dy*dy)/(2*parameter*parameter));
    }
    return mixture;
}

