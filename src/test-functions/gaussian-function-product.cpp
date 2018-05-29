#include "test-functions.h"

#include <math.h>
#include <iostream>

double gaussian_halfspace_product_points2d(double x, double y, double parameter, double angle,
                                           double xcenter, double ycenter,
                                     double minval, double maxval, double* domain){

    double gauss = gaussian_points2d(x, y, xcenter, ycenter, parameter);
    double halfspace = halfspace_points2d(x, y, angle, minval, maxval, domain);

    return gauss * halfspace;
}

double gaussian_disk_product_points2d(double x, double y, double parameter, double diskRadius,
                                      double centerx, double centery, double minvalue, double maxvalue){

    double gauss = gaussian_points2d(x, y, centerx, centery, parameter);
    double disk = disk_points2d(x, y, diskRadius, centerx, centery, minvalue, maxvalue);

    return gauss * disk;
}

double gaussian_box_product_points2d(double x, double y, double parameter, double asideLength, double bsideLength,
                                     double xcenter, double ycenter, double minvalue, double maxvalue,
                                     double boxAngle){

    double gauss = gaussian_points2d(x, y, xcenter, ycenter, parameter);
    double box = box_points2d(x, y, asideLength, bsideLength, minvalue, maxvalue, boxAngle);

    return gauss * box;
}

//double gaussian_gaussian_product_points2d(double x, double y, double sigma1, double sigma2){

//    double gauss1 = gaussian_points2d(x, y, sigma1);
//    double gauss2 = gaussian_points2d(x, y, sigma2);

//    return gauss1 * gauss2;
//}

