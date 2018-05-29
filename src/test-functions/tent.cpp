#include "test-functions.h"
#include "./../core/constants.h"
#include <math.h>
#include <iostream>

double tent_points2d(double x, double y, double param){

    //x = x + y * 1.4142;
   //0.8 - fabs(x);
    //return domain[2] - fabs(param * x);
    return 1.0 - fabs(param * x);
}

