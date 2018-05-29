#include "test-functions.h"
#include "./../core/constants.h"
#include <math.h>
#include <iostream>


double step_points1d(double x, double xhalflength, double minvalue, double maxvalue){

    if(x <= 0)
        return maxvalue;
    else
        return minvalue;
}

double box_points2d(double x, double y, double xhalflength, double yhalflength,
                    double minval, double maxval, double theta){

    if(theta == 0){
        if(x <= xhalflength && x >= -xhalflength && y <= yhalflength && y >= -yhalflength)
            return maxval;
        else
            return minval;
    }
    else if(theta == 90 || theta == 180 || theta == 270 || theta == 360){
        std::cerr << "Axis aligned boxes are already considered, use theta==0 !!!" << std::endl;
        exit(-1);
    }
    else{

        double angle = theta * deg2rad;
        double rotmat[] = {cos(angle), -sin(angle), sin(angle), cos(angle)};

        double xbox[] = {-xhalflength,-xhalflength,xhalflength,xhalflength};
        double ybox[] = {-yhalflength,yhalflength,yhalflength,-yhalflength};

        double xrotated[4], yrotated[4];
        for(int i=0; i < 4; i++){
            xrotated[i] = xbox[i] * rotmat[0] + ybox[i] * rotmat[1];
            yrotated[i] = xbox[i] * rotmat[2] + ybox[i] * rotmat[3];
        }

        double aline[] = {xrotated[0], yrotated[0], xrotated[1],yrotated[1]};
        double bline[] = {xrotated[1], yrotated[1], xrotated[2],yrotated[2]};
        double cline[] = {xrotated[2], yrotated[2], xrotated[3],yrotated[3]};
        double dline[] = {xrotated[3], yrotated[3], xrotated[0],yrotated[0]};

        double aslope = (aline[3] - aline[1]) / (aline[2] - aline[0]);
        double bslope = (bline[3] - bline[1]) / (bline[2] - bline[0]);
        double cslope = (cline[3] - cline[1]) / (cline[2] - cline[0]);
        double dslope = (dline[3] - dline[1]) / (dline[2] - dline[0]);

        double xyaline = (y - aline[1]) - (aslope * (x - aline[0]));
        double zeroaline = - aline[1] + (aslope * aline[0]);
        double asignproduct = xyaline * zeroaline;

        double xybline = (y - bline[1]) - (bslope * (x - bline[0]));
        double zerobline = - bline[1] + (bslope * bline[0]);
        double bsignproduct = xybline * zerobline;

        double xycline = (y - cline[1]) - (cslope * (x - cline[0]));
        double zerocline = - cline[1] + (cslope * cline[0]);
        double csignproduct = xycline * zerocline;

        double xydline = (y - dline[1]) - (dslope * (x - dline[0]));
        double zerodline = - dline[1] + (dslope * dline[0]);
        double dsignproduct = xydline * zerodline;

        ///Since (0,0) and (x,y) are on the same side of each line,
        /// return maxval only when all signproducts are positive.
        if( asignproduct >= 0 && bsignproduct >= 0 && csignproduct >= 0 && dsignproduct >= 0 ){
            return maxval;
        }
        else
            return minval;
    }
}


