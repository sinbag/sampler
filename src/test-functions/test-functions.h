#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H

#include "./../core/constants.h"
#include <iostream>

double step_points1d(double x, double xhalflength, double minvalue=0.0, double maxvalue=1.0);
double disk_points2d(double x, double y, double diskRadius, double centerx, double centery,
                     double minvalue=0.0, double maxvalue=1.0);

double box_points2d(double x, double y, double xhalflength, double yhalflength,
                    double minval, double maxval, double theta=0);

double tent_points2d(double x, double y, double param);

double constant_linear_constant(double x, double y, double* domain);
double constant_linear(double x, double y, double* domain);
double linear(double x, double y, double* domain);
double plane(double x, double y, double* domain);
double linear_linear(double x, double y, double* domain);
double triangular_distribution(double x, double xPeakPos, double* domain);

double rotated_constant_linear_constant(double x, double y, double* domain);
double quadratic(double x, double y, double* domain);

double c2_discontinuity(double x, double y, double param, double* domain);
double c1_discontinuity(double x, double y, double param, double* domain);

double gaussian_points2d(double x, double y, double cx, double cy, double parameter);
double gaussian_points1d(double x, double xcenter, double parameter);

double disk_gaussian_points2d(double x, double y, double x1, double y1, double xcenter, double ycenter, double parameter);
double disk_gaussian_points1d(double x, double x1, double xcenter, double parameter);

double halfspace_points2d(double ptx, double pty, double halfSpaceAngle, double minval, double maxval, double* domain);
double halfspace(double x, double y, double param, double maxvalue, double minvalue, double* domain);
double gaussian_halfspace_product_points2d(double x, double y, double parameter, double angle,
                                           double xcenter, double ycenter,
                                           double minvalue, double maxvalue, double *domain);
double gaussian_disk_product_points2d(double x, double y, double parameter, double diskRadius,
                                      double centerx, double centery, double minvalue, double maxvalue);
double gaussian_box_product_points2d(double x, double y, double parameter, double asideLength,
                                     double bsideLength, double xcenter, double ycenter,
                                     double minvalue, double maxvalue, double boxAngle=0.0);

double gaussian_gaussian_product_points2d(double x, double y, double sigma1, double sigma2);

class ck_continuity_irrational{
public:
    double constant_linear_constant(double x, double y, double* domain);
    double constant_linear(double x, double y, double* domain);
    double linear(double x, double y, double* domain);
    double linear_linear(double x, double y, double* domain);
};

#endif

