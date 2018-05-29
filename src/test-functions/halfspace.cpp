#include "test-functions.h"
#include <math.h>
#include <iostream>

double halfspace_points2d(double ptx, double pty, double halfSpaceAngle, double minval, double maxval, double* domain){

    ///line coordinates of the boundary of the halfspace
    double lxi = 0.0, lyi = 0.0, lxf = 0.0, lyf = 0.0;

    ///We consider only lines with slopes in this range
    if(halfSpaceAngle >=0.0 && halfSpaceAngle <= 45.0){
        lxi = domain[0];
        lyi = -tan(halfSpaceAngle * deg2rad) * fabs(domain[2]);
        lxf = domain[2];
        lyf = tan(halfSpaceAngle * deg2rad) * fabs(domain[2]);
    }
    else{
        std::cerr << "Please enter angle between 0 and 45 for halfspace !!!" << std::endl;
        exit(-1);
    }

    double DX = lxf - lxi;
    double DY = lyf - lyi;

    ///Check which end-point of the input line belongs to the non-zero halfspace.
    double flag =0.0;
    if(DX == 0){
        /// y - slope * x = 0;
        flag = ptx-lxi;
    }
    else{
        double m = DY / DX;
        flag =  (pty-lyi) - (m * (ptx-lxi)); ///first end-point of the input line
    }

    if(flag > 0)
        return maxval;
    else
        return minval;
}

double halfspace(double x, double y, double param, double maxvalue, double minvalue, double* domain){

    double flag = x - 1.414*y + param;
    if (flag < 0)
        return maxvalue;
    else
        return minvalue;
}

