#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <random>

#include "./core/utils.h"
#include "./io/write-exr.h"
#include "./io/write-eps.h"
#include "core/constants.h"
#include "test-functions/test-functions.h"
#include "statistics/statistics.h"
#include "util-samples/util-samples.h"

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
    unsigned int numTrials, startIndex, endIndex;
    double parameter  = 0.0;
    std::string testfunction = "";
    std::string mode;
    std::string wrapper;
    std::string samplingpattern = "bnot";

    unsigned short int seed;

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
            ("numTrials,n",
             boostPO::value<unsigned int>(&numTrials)->default_value(200),
             "Number of trials")
            ("initIndex,i",
             boostPO::value<unsigned int>(&startIndex)->default_value(1),
             "start from given number of samples")
            ("endIndex,e",
             boostPO::value<unsigned int>(&endIndex)->default_value(1024),
             "end at given number of samples")
            ("function,f",
             boostPO::value<std::string>(&testfunction)->default_value("disk"),
             "Test function")
            ("parameter,p",
             boostPO::value<double>(&parameter)->default_value(0.25),
             "Function parameter")
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

    //============================================================
    Sampler sampler(fn_rules, fn_bary, fn_dev);
    //============================================================

    std::stringstream ss;
    //============================================================
    boost::filesystem::path source_dir_path( boost::filesystem::current_path() );

    std::string datafiles, images, graphs;
    ss.str(std::string());
    ss << source_dir_path.string() << "/results/";
    std::string resultFolder = ss.str();
    mkdir(resultFolder.c_str() ,0755);
    ss.str(std::string());
    ss << resultFolder << "variance-analysis-points-polyhex-bnot" << "/";

    mkdir(ss.str().c_str(),0755);

    create_folders(ss.str(), datafiles, images, graphs);
    //============================================================

    double bBoxMin = -0.5;
    double bBoxMax = 0.5;
    double domain[] = {bBoxMin, bBoxMin, bBoxMax, bBoxMax};

    //##########################################################
    double cx=0;
    double cy=0;
    std::ofstream fvar, fmean;
    ss.str(std::string());
    ss << datafiles << "variance-pointsamples-" << testfunction << "-cx" << cx << "-cy" << cy << "-bbmin" << bBoxMin << "-bbmax" << bBoxMax << "-" << mode << "-polyhex" << samplingpattern  << ".txt";
    fvar.open(ss.str().c_str(), std::ios::app);

    ss.str(std::string());
    ss << datafiles << "mean-pointsamples-" << testfunction  << "-cx" << cx << "-cy" << cy << "-bbmin" << bBoxMin << "-bbmax" << bBoxMax << "-" << mode << "-polyhex" << samplingpattern  << ".txt";
    fmean.open(ss.str().c_str(), std::ios::app);
    //##########################################################

    std::random_device rd;
    static thread_local std::mt19937 generator(rd());
    std::uniform_int_distribution<int> dis(0, std::numeric_limits<int>::max());
    int randomseed = dis(generator);
    
    if( seed == 0 )
    {
        srand48(time(NULL));
    }
    if( vm.count("seed") ) seed = (seed-1)%408;

    for(int N = startIndex; N <= endIndex;){

        int totalSamples = N*N;

        double mean = 0.0;
        double variance = 0.0;

    //WriterFileRaw write(fn_out);
    for(int trial = 1; trial <= numTrials; trial++){

        seed = std::ceil(drand48()*408);
        WriterFileRaw write(fn_out);
//        WriterVector write;
        sampler.generateUniform(totalSamples, 0, write, seed);

        //##########################################################

        totalSamples = write.pts().size() * 0.5;
        fprintf(stderr,"\r N trial (%d %d)", totalSamples, trial);

        //##########################################################

        std::vector<float> pointset;

        for(int i =0; i < write.pts().size();i+=2){

            float x = write.pts().at(i);//.x();
            float y = write.pts().at(i+1);//.y();

            pointset.push_back(x + domain[0]);
            pointset.push_back(y + domain[1]);
            //std::cout << x << " "<< y << std::endl;
        }

        //##########################################################

        std::vector<float> finalsamples;
        if(mode == "homogenize")
            finalsamples = homogenize_pointsamples(pointset, domain, 2);
        else
            finalsamples = pointset;
        //##########################################################

        totalSamples = finalsamples.size() * 0.5;

        double integral = 0;
        for(int i = 0; i < totalSamples; i++){
            double x = finalsamples[2*i+0];
            double y = finalsamples[2*i+1];
            
            if(testfunction == "gaussian")
                integral += gaussian_points2d(x,y, 0.1,0.1,parameter);
            else if(testfunction == "disk")
                integral += disk_points2d(x, y, parameter, 0., 0.);
            else if(testfunction == "tent")
                integral += tent_points2d(x,y, parameter);
            else{
                std::cerr << "WARNING: Integrand not declared !!!" << std::endl;
            }
        }

        double measure = (domain[2]-domain[0])*(domain[3]-domain[1]);
        integral *= measure;
        integral /= double(totalSamples);

        //##########################################################

        //Incremental Mean and Variance computation
        mean = mean_progressive(mean, integral, trial);
        variance = variance_progressive(variance, mean, integral, trial);
    }

    fmean << std::fixed << totalSamples << std::setprecision(20) << " " << mean << std::endl;
    fvar << std::fixed << totalSamples << std::setprecision(20) << " " << variance << std::endl;

    if(N < 31)
        N++;
    else
        N *=2;
    }
    fvar.close();
    fmean.close();
    std::cerr << std::endl;
    return 0;
}

