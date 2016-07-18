#include "fourier-analysis.h"
#include <iostream>
#include <fstream>
#include <iomanip>

template <typename T>
void FT<T>::realimag_fourier_spectrum(T *realSpectrum, T* imagSpectrum,
                                         const std::complex<T>* complexSpectrum,
                                         int width, int height){
    for(int r = 0; r < height; r++){
        for(int c = 0; c < width; c++){
            int index = r*width+c;
            T real = complexSpectrum[(r*width+c)].real();
            T imag = complexSpectrum[(r*width+c)].imag();

            realSpectrum[index] = real;
            imagSpectrum[index] = imag;
        }
    }
}

template <typename T>
void FT<T>::power_fourier_spectrum(T *powerSpectrum, const std::complex<T>* complexSpectrum, int totalSamples, int width, int height){
    for(int r = 0; r < height; r++){
        for(int c = 0; c < width; c++){
            T real = complexSpectrum[(r*width+c)].real();
            T imag = complexSpectrum[(r*width+c)].imag();

            T power = (real*real + imag*imag) / (totalSamples);

            powerSpectrum[(r*width+c)] = power;
        }
    }
}

template <typename T>
void FT<T>::magnitude_fourier_spectrum(T *magnitudeSpectrum, const std::complex<T>* complexSpectrum, int totalSamples, int width, int height){
    for(int r = 0; r < height; r++){
        for(int c = 0; c < width; c++){
            T real = complexSpectrum[(r*width+c)].real();
            T imag = complexSpectrum[(r*width+c)].imag();

            T power = (real*real + imag*imag) / (totalSamples);

            magnitudeSpectrum[(r*width+c)] = sqrt(power);
        }
    }
}

template <typename T>
void FT<T>::phase_fourier_spectrum(T* phaseSpectrum, const std::complex<T> *complexSpectrum, int width, int height){

    for(int r = 0; r < height; r++){
        for(int c = 0; c < width; c++){
            T real = complexSpectrum[(r*width+c)].real();
            T imag = complexSpectrum[(r*width+c)].imag();

            T phase = atan2(imag, real);

            phaseSpectrum[(r*width+c)] = phase;
        }
    }
}

template <typename T>
void FT<T>::compute_radial_mean_powerspectrum(std::string filename, T *power,
                                              int width, int height){

    if(width != height){
        std::cerr << "We assume square images for radial mean power computation !!!" << std::endl;
        exit(-2);
    }
    ///Radial Power spectrum
    int halfwidth = width*0.5;

    T* radialHistogram = new T[halfwidth]();
    int* histoCounter = new int[halfwidth]();
    int xcenter = halfwidth;
    int ycenter = halfwidth;
    for(int r = 0; r < height; r++){
        for(int c = 0; c < width; c++){
            double dx = xcenter-c;
            double dy = ycenter-r;
            double distance = sqrt(dx*dx+dy*dy);

            int imgIndex = r*width+c;
            int index = distance;
            if(distance > halfwidth-1)
                continue;
            else{
                radialHistogram[index] += power[imgIndex];
                histoCounter[index] += 1;
            }
        }
    }

    for(int i = 0; i < halfwidth; i++)
        radialHistogram[i] /= double(histoCounter[i]);

    std::ofstream file;
    file.open(filename);

    for(int i = 0; i < halfwidth-5; i++)
        file << i << " " << std::fixed << std::setprecision(15) <<  radialHistogram[i] << std::endl;
    file.close();


    delete [] histoCounter;
    delete [] radialHistogram;
}


template void FT<double>::magnitude_fourier_spectrum(double* magnitudeSpectrum, const std::complex<double>* complexSpectrum, int totalSamples, int width, int height);
template void FT<double>::power_fourier_spectrum(double* powerSpectrum, const std::complex<double>* complexSpectrum, int totalSamples, int width, int height);
template void FT<double>::phase_fourier_spectrum(double* phaseSpectrum, const std::complex<double>* complexSpectrum, int width, int height);
template void FT<double>::realimag_fourier_spectrum(double* realSpectrum, double* imagSpectrum, const std::complex<double>* complexSpectrum, int width, int height);


template void FT<float>::magnitude_fourier_spectrum(float* magnitudeSpectrum, const std::complex<float>* complexSpectrum, int totalSamples, int width, int height);
template void FT<float>::power_fourier_spectrum(float* powerSpectrum, const std::complex<float>* complexSpectrum, int totalSamples, int width, int height);
template void FT<float>::phase_fourier_spectrum(float* phaseSpectrum, const std::complex<float>* complexSpectrum, int width, int height);

template void FT<double>::compute_radial_mean_powerspectrum(std::string filename, double *power, int width, int height);
template void FT<float>::compute_radial_mean_powerspectrum(std::string filename, float *power, int width, int height);
