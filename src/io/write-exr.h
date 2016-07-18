#ifndef WRITE_EXR_H
#define WRITE_EXR_H
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <ImfRgbaFile.h>
#include <half.h>
#include <assert.h>

template <typename T>
void write_exr_grey(std::string name, T* pixels, int xRes, int yRes);
/*
template <typename T>
void write_exr_float(std::string name, T pixels, int xRes, int yRes);
template <typename T>
void write_exr_rgba(std::string name, T pixels, int xRes, int yRes);
template <typename T>
void write_exr_rgb(std::string name, T pixels, int xRes, int yRes);
*/
#endif
