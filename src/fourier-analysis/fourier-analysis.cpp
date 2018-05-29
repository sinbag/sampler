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
void FT<T>::power_fourier_spectrum(T *powerSpectrum, const std::complex<T>* complexSpectrum, int totalSamples, int width, int height, int ndims){
    if(ndims == 2){
        for(int r = 0; r < height; r++){
            for(int c = 0; c < width; c++){
                T real = complexSpectrum[(r*width+c)].real();
                T imag = complexSpectrum[(r*width+c)].imag();
                
                T power = (real*real + imag*imag) / (totalSamples);
                
                powerSpectrum[(r*width+c)] = power;
            }
        }
    }
    else if(ndims == 1){
        for(int c = 0; c < width; c++){
            T real = complexSpectrum[(c)].real();
            T imag = complexSpectrum[(c)].imag();
            
            T power = (real*real + imag*imag) / (totalSamples);
            
            powerSpectrum[(c)] = power;
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
void FT<T>::radialStatistics(std::string meanfilename,
                              std::string anisotropyfilename,
                              T *power,int width, int height,
                              int trialCount, int ndims){
                                              

//    if(width != height){
//        std::cerr << "We assume square images for radial mean power computation !!!" << std::endl;
//        exit(-2);
//    }
    ///Radial Power spectrum
    int halfwidth = width*0.5;

    std::vector<T> radialHistogram(halfwidth, T(0.0));
    std::vector<T> radialVariance(halfwidth, T(0.0));
    std::vector<T> radialAnisotropy(halfwidth, T(0.0));
//    T* radialVariance = new T[halfwidth]();
//    T* radialAnisotropy = new T[halfwidth]();
    int* histoCounter = new int[halfwidth]();
    
    int xcenter = halfwidth;
    int ycenter = halfwidth;
    if(ndims == 2)
    {
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
    }
    else if(ndims == 1){
        for(int r = xcenter; r < width; r++){
            int index = r - xcenter;
            {
                //std::cerr << r << " " << index << std::endl;
                radialHistogram[index] += power[r];
                histoCounter[index] += 1;
            }
        }
    }

    for(int i = 0; i < halfwidth; i++){
        radialHistogram[i] /= double(histoCounter[i]);
        histoCounter[i] = 0.0;
    }

    ///Anisotropy
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
                double deviation = power[imgIndex] - radialHistogram[index];
                radialVariance[index] += (deviation*deviation);
                histoCounter[index] += 1;
            }
        }
    }
    
    for(int k = 1; k < halfwidth; k++){
        radialVariance[k] /= (histoCounter[k]-1);
        radialAnisotropy[k] = 10*log10(radialVariance[k] / (radialHistogram[k]*radialHistogram[k]));
        
        ///
        /// Anisotropy scales by log10(numTrials), therefore, we divide by this
        /// number below to make sure Anisotropy always scales to -10dB.
        ///
        if(trialCount > 1)
            radialAnisotropy[k] /= log10(trialCount);
    }
    
    std::ofstream mfile, afile;
    mfile.open(meanfilename.c_str());
    afile.open(anisotropyfilename.c_str());
    
    for(int i = 0; i < halfwidth; i++){
        mfile << i << " " << std::fixed << std::setprecision(15) <<  radialHistogram[i] << std::endl;
        afile << i << " " << std::fixed << std::setprecision(15) <<  radialAnisotropy[i] << std::endl;
    }
    mfile.close();
    afile.close();
}


template void FT<double>::magnitude_fourier_spectrum(double* magnitudeSpectrum, const std::complex<double>* complexSpectrum, int totalSamples, int width, int height);
template void FT<double>::power_fourier_spectrum(double* powerSpectrum, const std::complex<double>* complexSpectrum, int totalSamples, int width, int height,int ndims);
template void FT<double>::phase_fourier_spectrum(double* phaseSpectrum, const std::complex<double>* complexSpectrum, int width, int height);
template void FT<double>::realimag_fourier_spectrum(double* realSpectrum, double* imagSpectrum, const std::complex<double>* complexSpectrum, int width, int height);


template void FT<float>::magnitude_fourier_spectrum(float* magnitudeSpectrum, const std::complex<float>* complexSpectrum, int totalSamples, int width, int height);
template void FT<float>::power_fourier_spectrum(float* powerSpectrum, const std::complex<float>* complexSpectrum, int totalSamples, int width, int height, int ndims);
template void FT<float>::phase_fourier_spectrum(float* phaseSpectrum, const std::complex<float>* complexSpectrum, int width, int height);

template void FT<double>::radialStatistics(std::string mfilename,std::string afilename,
                                           double *power, int width, int height, int trial, int ndims);
template void FT<float>::radialStatistics(std::string mfilename,std::string afilename,
                                           float *power, int width, int height, int trial, int ndims);
