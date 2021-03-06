#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <sys/stat.h>

#include "./core/utils.h"
#include "./io/write-exr.h"
#include "./io/write-eps.h"
#include "./fourier-analysis/fourier-analysis.h"

#include "sampler.hpp"
#include "tiling.hpp"

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
    float* power = new float[width*height]();
    float* powerAccum = new float[width*height]();
    std::complex<float>* complexSpectrum = new std::complex<float>[width*height]();
    int N = nbSample;
    float spaceSize = 0.21;
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
//        std::cerr << "# of samples: " << write.pts().size() << " "<< nbSample << std::endl;
        //##########################################################

        for(int k = 0; k < width*height; k++){
            power[k] = 0;
        }
        //##########################################################

        FT<float>::continuous_fourier_spectrum_parallel(complexSpectrum, finalsamples, width, height, frequencyStep, ndims, projection);
        FT<float>::power_fourier_spectrum(power, complexSpectrum, nbSample, width, height, ndims);
        //perform_cft_parallel(power, finalsamples, N, width, height,frequencyStep);

        //##########################################################

        for(int i = 0; i < width*height;i++)
            powerAccum[i] += power[i];

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
            ss << images << "power-cft-fstep" << frequencyStep << "-" << mode << "-" << samplingpattern << "-" << wrapper << "-n" << N << "-r" << resolution << "-" << s1 << ".png";
            write_exr_grey(ss.str(), power, width, height);

            ss.str(std::string());
            ss << images << "pointset-" << mode << "-" << samplingpattern<< "-" << wrapper << "-n" << N << "-r" << resolution << "-" << s1 << ".png";
            write_eps(ss.str(), finalsamples);
            write_txt(ss.str(), finalsamples);
            //##########################################################
            
            ss.str(std::string());
            ss << datafiles << "radial-mean-cft-fstep" << frequencyStep << "-"  << mode << "-" << samplingpattern << "-" << wrapper << "-n" << N << "-r" << resolution << "-proj" << projection << "-" << s1 << ".txt";
            std::string meanfilename = ss.str();
            
            ss.str(std::string());
            ss << datafiles << "radial-anisotropy-cft-fstep" << frequencyStep << "-"  << mode << "-" << samplingpattern << "-" << wrapper << "-n" << N << "-r" << resolution << "-proj" << projection << "-" << s1 << ".txt";
            std::string anisotropyfilename = ss.str();
            
            FT<float>::radialStatistics(meanfilename, anisotropyfilename, power,
                                        width, height, patch, ndims);
        }
    }
    std::cerr << std::endl;

    delete [] power;
    delete [] powerAccum;
    delete [] complexSpectrum;
}


