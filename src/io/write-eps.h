#ifndef WRITEEPS_H
#define WRITEEPS_H

#include <iostream>
#include <vector>

template <typename T>
void write_eps(std::string filename, const std::vector<T>& points, int nDims=2,
               double radius=2.0, double scale=512.0);

template <typename T>
void write_txt(std::string filename, const std::vector<T>& points, int nDims=2);

#endif // WRITEEPS_H
