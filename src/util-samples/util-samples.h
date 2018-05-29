#ifndef UTILSAMPLES_H
#define UTILSAMPLES_H

#include <iostream>
#include <vector>

std::vector<double> homogenize_pointsamples2d(std::vector<double> &vec, double* domain);

template <typename T>
std::vector<T> homogenize_pointsamples(std::vector<T> &vec, double* domain, int nDims);

#endif // UTILSAMPLES_H
