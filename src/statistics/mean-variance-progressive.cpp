#include "statistics.h"
#include <iostream>

double mean_progressive(double &mean, const double &integral, int trial){

    if(trial == 0){
        std::cerr <<"Progressive mean trial index must start from 1, and not 0 !!!" << std::endl;
        exit(-2);
    }
    else{
        mean = ((trial-1)/double(trial)) * mean + (1/double(trial)) * integral;
    }

    return mean;
}

double variance_progressive(double &variance, const double &mean, const double &integral, int trial){

    if(trial == 0){
        std::cerr <<"Progressive variance trial index must start from 1, and not 0 !!!" << std::endl;
        exit(-2);
    }
    else if(trial < 2){
        variance = 0;
    }
    else{
    variance = ((trial-1)/double(trial)) * variance + (1/double(trial-1)) * (integral - mean) * (integral - mean);
    }
    return variance;
}

//double variance_online
