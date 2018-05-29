//
//  cft-cuda.cpp
//  fourier-analysis
//
//  Created by Gurprit Singh on 5/29/18.
//

#include <iostream>

#ifdef USE_CUDA_ENABLED
#include <utils_cuda.h>
#include <timer_cuda.h>


#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <sys/stat.h>

#include "./core/constants.h"
#include "./core/utils.h"
#include "./io/write-exr.h"
#include "./io/write-eps.h"
#include "./fourier-analysis/fourier-analysis.h"

#include "sampler.hpp"
#include "tiling.hpp"

void cftmd_cuda(double * d_outReal, double * d_outImag,
                double* d_inPoints,  double* d_inFuncVals,
                int npts, int resolution,
                int ndims, double twoPiConst, float frequencyStep);


namespace boostPO = boost::program_options;

int main(int argc, char** argv)
{
    /* ARG PARSER *****************************************************/
    std::string fn_rules;
    std::string fn_bary;
    std::string fn_dev;
    std::string fn_out="test.txt";
    unsigned int nbSample, numPatches, trialStep,resolution, ndims, projection;
    float frequencyStep;
    std::string mode;
    std::string wrapper;
    std::string samplingpattern = "bnot";
    
    unsigned short int seed;
    std::cerr <<
    "usage: ./fourier-analysis -r ./../data/lut/production_rules.dat -b ./../data/lut/barycenters.dat -d ./../data/lut/offsets_bnot.dat -n 1024 -o sampling.txt" << std::endl;
    boostPO::variables_map vm;
    boostPO::options_description desc("Allowed options");
    desc.add_options()
    ("help,h",
     "produce help message")
    ("rule,r",
     boostPO::value<std::string>(&fn_rules)->required(),
     "REQUIRED | Subdivision rule filename")
    ("bary,b",
     boostPO::value<std::string>(&fn_bary)->required(),
     "REQUIRED | Barycenter offset filename (for each rule id)")
    ("dev,d",
     boostPO::value<std::string>(&fn_dev),
     "Offset LUT filename (for each structural indices)")
    ("nbSample,n",
     boostPO::value<unsigned int>(&nbSample)->default_value(1024),
     "Number of sample de generate")
    ("patches,p",
     boostPO::value<unsigned int>(&numPatches)->default_value(200),
     "Number of patches")
    ("avg,a",
     boostPO::value<unsigned int>(&trialStep)->default_value(1),
     "Output average after trialStep patches")
    ("resolution,u",
     boostPO::value<unsigned int>(&resolution)->default_value(512),
     "resolution")
    ("freq,f",
     boostPO::value<float>(&frequencyStep)->default_value(1.0),
     "Frequency Step")
    ("dims,e",
     boostPO::value<unsigned int>(&ndims)->default_value(2),
     "Number of dimensions")
    ("projection,c",
     boostPO::value<unsigned int>(&projection)->default_value(0),
     "projection")
    ("translate,t",
     boostPO::value<std::string>(&mode)->default_value("nohomogenize"),
     "Homogenization mode")
    ("wrapper,w",
     boostPO::value<std::string>(&wrapper)->default_value("toroidal"),
     "Toroidal or nontoroidal")
    ("seed,s",
     boostPO::value<unsigned short int>(&seed)->default_value(0),
     "Initial tile to use for sampling ([1-408], 0 = random)")
    ;
    
    try
    {
        boostPO::store(
                       boostPO::command_line_parser(argc, argv).
                       options(desc).run(), vm);
        boostPO::notify(vm);
    }
    catch(boost::program_options::error& e)
    {
        std::cerr << e.what() << std::endl;
        std::cout << desc << "\n";
        exit(EXIT_FAILURE);
    }
    
    if(vm.count("help"))
    {
        std::cout << desc << "\n";
        exit(EXIT_SUCCESS);
    }
    
    /* PROG ***********************************************************/
    Sampler sampler(fn_rules, fn_bary, fn_dev);
    
    std::stringstream ss;
    //============================================================
    boost::filesystem::path source_dir_path( boost::filesystem::current_path() );
    
    std::string datafiles, images, graphs;
    ss.str(std::string());
    ss << source_dir_path.string() << "/results/";
    std::string resultFolder = ss.str();
    mkdir(resultFolder.c_str() ,0755);
    ss.str(std::string());
    ss << resultFolder << "powerspectrum-cft-bnot-n" << nbSample << "/";
    
    mkdir(ss.str().c_str(),0755);
    
    create_folders(ss.str(), datafiles, images, graphs);
    //============================================================
    
    int width = resolution, height = resolution;
    if(ndims == 1)
        height = 1;
    int _numPixels = width*height;
    int _ndims = 2;
    float* power = new float[width*height]();
    float* powerAccum = new float[width*height]();
    std::complex<float>* complexSpectrum = new std::complex<float>[width*height]();
    int tempN = nbSample;
    float spaceSize = 0.21;
    //============================================================
    
    double* h_outReal = new double[_numPixels]();
    double* h_outImag = new double[_numPixels]();
    
    //allocate memory on the device for both input and output
    checkCudaErrors(cudaMalloc(&d_outReal, sizeof(double) * _numPixels));
    checkCudaErrors(cudaMalloc(&d_outImag, sizeof(double) * _numPixels));
    
    //Initializing device storage here
    checkCudaErrors(cudaMemset(d_outReal, 0, sizeof(double) * _numPixels));
    checkCudaErrors(cudaMemset(d_outImag, 0, sizeof(double) * _numPixels));
    
    //============================================================
    
    if( seed == 0 )
    {
        srand48(time(NULL));
    }
    if( vm.count("seed") ) seed = (seed-1)%408;
    
    //WriterFileRaw write(fn_out);
    for(int patch = 1; patch <= numPatches; patch++){
        
        seed = std::ceil(drand48()*408);
        WriterFileRaw write(fn_out);
        //        WriterVector write;
        sampler.generateUniform(N, 0, write, seed);
        
        //##########################################################
        
        nbSample = write.pts().size() * 0.5;
        fprintf(stderr,"\r N trial (%d %d)", nbSample, patch);
        
        //##########################################################
        
        std::vector<float> finalsamples;
        
        for(int i =0; i < write.pts().size();i+=2){
            
            float x = write.pts().at(i);//.x();
            float y = write.pts().at(i+1);//.y();
            finalsamples.push_back(x);
            finalsamples.push_back(y);
            //            std::cout << x << " "<< y << std::endl;
        }
        
        int N = finalsamples.size() * 0.5;
        //============================================================
        ///Host CPU array to store points to pass to GPU via CUDA
        double* h_inPointsArray = new double[_ndims*N]();
        double* h_inFuncVals = new double[N]();
        
        //allocate memory on the device for both input and output
        checkCudaErrors(cudaMalloc(&d_inPoints, _ndims * sizeof(double) * N));
        checkCudaErrors(cudaMalloc(&d_inFuncVals, sizeof(double) * N));
        
        checkCudaErrors(cudaMemset(d_inPoints, 0, _ndims * sizeof(double) * N));
        checkCudaErrors(cudaMemset(d_inFuncVals, 0, sizeof(double) * N));
        
        for(int i=0; i < tempN; i++){
            h_inPointsArray[i] = x;
            h_inPointsArray[i+1] = y;
            h_inFuncVals[i] = 1;
        }

        //============================================================
        
        //        std::cerr << "# of samples: " << write.pts().size() << " "<< nbSample << std::endl;
        //##########################################################
        
        for(int k = 0; k < width*height; k++){
            power[k] = 0;
            h_outReal[i] = 0.;
            h_outImag[i] = 0.;
        }
        
        //============================================================
        ///Copying sample locations to device
        checkCudaErrors(cudaMemcpy(d_inPoints, h_inPointsArray, _ndims * sizeof(double) * N, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(d_inFuncVals, h_inFuncVals, sizeof(double) * N, cudaMemcpyHostToDevice));
        
#ifdef GPU_TIMER
        GpuTimer timer;
        timer.Start();
        //            continuous_fourier_transform_cuda(d_outReal, d_outImag, d_inPoints,
        //                                              d_inFuncVals, n, _xRes, _xRes,
        //                                              twopi, _frequencyStep);
        cftmd_cuda(d_outReal, d_outImag, d_inPoints,
                   d_inFuncVals, N, width, _ndims,
                   twopi, 1.0);
        timer.Stop();
        cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
        
        int err = printf("\nYour code ran in: %f msecs.\n", timer.Elapsed());
        if (err < 0) {
            //Couldn't print! Probably the student closed stdout - bad news
            std::cerr << "Couldn't print timing information! STDOUT Closed!" << std::endl;
            exit(1);
        }
#else //GPU_TIMER
        //            continuous_fourier_transform_cuda(d_outReal, d_outImag, d_inPoints,
        //                                              d_inFuncVals, n, _xRes, _xRes,
        //                                              twopi, _frequencyStep);
        cftmd_cuda(d_outReal, d_outImag, d_inPoints,
                   d_inFuncVals, N, width, _ndims,
                   twopi, 1.0);
        cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
#endif //ENDIF GPU_TIMER
        
        ///Copy computed Fourier coeffs from GPU back to host arrays
        checkCudaErrors(cudaMemcpy(h_outReal, d_outReal, sizeof(double) * _numPixels, cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(h_outImag, d_outImag, sizeof(double) * _numPixels, cudaMemcpyDeviceToHost));
        
        ///Copy Fourier coeffs to _complexSpectrum array
        for(int i=0; i  < _numPixels; i++){
            complexSpectrum[i].real(h_outReal[i]);
            complexSpectrum[i].imag(h_outImag[i]);
        }
        //============================================================
        
        
        //##########################################################
        
        FT<float>::power_fourier_spectrum(power, complexSpectrum, nbSample, width, height, ndims);
        //perform_cft_parallel(power, finalsamples, N, width, height,frequencyStep);
        
        //##########################################################
        
        for(int i = 0; i < width*height;i++){
            powerAccum[i] += power[i];
        }
        
        if(patch % trialStep == 0 || patch == 1){
            
            ss.str(std::string());
            ss << patch;
            std::string s1 = ss.str();
            paddedzerosN(s1, numPatches);
            //##########################################################
            
            for(int i = 0; i < width*height; i++)
                power[i] = powerAccum[i] / patch;
            
            //##########################################################
            
            ss.str(std::string());
            ss << images << "power-cft-fstep" << frequencyStep << "-" << mode << "-" << samplingpattern << "-" << wrapper << "-n" << tempN << "-r" << resolution << "-" << s1 << ".png";
            write_exr_grey(ss.str(), power, width, height);
            
            ss.str(std::string());
            ss << images << "pointset-" << mode << "-" << samplingpattern<< "-" << wrapper << "-n" << tempN << "-r" << resolution << "-" << s1 << ".png";
            write_eps(ss.str(), finalsamples);
            write_txt(ss.str(), finalsamples);
            //##########################################################
            
            ss.str(std::string());
            ss << datafiles << "radial-mean-cft-fstep" << frequencyStep << "-"  << mode << "-" << samplingpattern << "-" << wrapper << "-n" << tempN << "-r" << resolution << "-proj" << projection << "-" << s1 << ".txt";
            std::string meanfilename = ss.str();
            
            ss.str(std::string());
            ss << datafiles << "radial-anisotropy-cft-fstep" << frequencyStep << "-"  << mode << "-" << samplingpattern << "-" << wrapper << "-n" << tempN << "-r" << resolution << "-proj" << projection << "-" << s1 << ".txt";
            std::string anisotropyfilename = ss.str();
            
            FT<float>::radialStatistics(meanfilename, anisotropyfilename, power,
                                        width, height, patch, ndims);
        }
        delete [] h_inPointsArray;
        delete [] h_inFuncVals;
        cudaFree(d_inPoints);
    }
    
    std::cerr << std::endl;
    
    delete [] power;
    delete [] powerAccum;
    delete [] complexSpectrum;
    
    delete [] h_outReal;
    delete [] h_outImag;
    
    cudaFree(d_outReal);
    cudaFree(d_outImag);

}
#endif
