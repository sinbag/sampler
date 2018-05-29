#ifndef SPECTRALANALYSIS_H
#define SPECTRALANALYSIS_H

#include <complex>
#include <vector>

template <typename T>
class FT{
public:

    static void continuous_fourier_spectrum_parallel(std::complex<T> *complexSpectrum, std::vector<T> &points, int width, int height, T dstep, int ndims=2, int projection=0);
    static void discrete_fourier_spectrum(std::complex<T> *complexSpectrum, const T *mydata, int width, int height);

    static void power_fourier_spectrum(T *powerSpectrum, const std::complex<T>*complexSpectrum, int total_samples, int width, int height,int ndims=2);

    static void magnitude_fourier_spectrum(T *magnitudeSpectrum, const std::complex<T>*complexSpectrum, int total_samples, int width, int height);

    static void phase_fourier_spectrum(T *phaseSpectrum, const std::complex<T>*complexSpectrum, int width, int height);

    ///Inverse using the complex spectrum
    static void inverse_discrete_fourier_transform(T *spatialSignal, const std::complex<T> *complexSpectrum, int width, int height);

    ///Inverse using the magnitude and the phase spectrum
    static void inverse_discrete_fourier_transform(T *spatialSignal, T *magnitudeSpectrum, T *phaseSpectrum, int totalSamples, int width, int height);

    static void realimag_fourier_spectrum(T *realSpectrum, T* imagSpectrum,
                                             const std::complex<T>* complexSpectrum, int width, int height);

    ///Radial / Anisotropy / Differential Tools
    static void radialStatistics(std::string meanfilename,
                                 std::string anisotropyfilename,
                                 T *power,int width, int height,
                                 int trialCount, int ndims=2);
//                                                  std::string filename, T *indata, int width, int height);

    static void swapQuadrants(T *img, int width, int height);
};

#ifdef USE_CUDA_ENABLED
double *d_outReal;
double *d_outImag;
double *d_inPoints;
double *d_inFuncVals;
#endif

#endif // SPECTRALANALYSIS_H
