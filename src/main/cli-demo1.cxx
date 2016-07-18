#include <boost/program_options.hpp>
#include <iostream>

#include "sampler.hpp"
#include "tiling.hpp"

namespace boostPO = boost::program_options;

int main(int argc, char** argv)
{
    /* ARG PARSER *****************************************************/
    std::string fn_rules;
    std::string fn_bary;
    std::string fn_dev;
    std::string fn_out;
    unsigned int nbSample;
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
        ("out,o",
         boostPO::value<std::string>(&fn_out),
         "Output filename")
        ("nbSample,n",
         boostPO::value<unsigned int>(&nbSample)->default_value(1024),
         "Number of sample de generate")
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

    if( seed == 0 )
    {
        srand48(time(NULL));
    }
    if( vm.count("seed") ) seed = (seed-1)%408;

    if( vm.count("out") )
    {
        //WriterFileRaw write(fn_out);
        for(int p = 0; p < 2; p++){
            seed = std::ceil(drand48()*408);
            WriterVector write;
            sampler.generateUniform(nbSample, -1, write, seed);

            for(int i =0; i < write.pts().size();i++){
                std::cout << write.pts().at(i).x() <<" " << write.pts().at(i).y() << std::endl;
            }
            std::cout << std::endl;
            std::cout << "New Patch " << std::endl;
            std::cout << std::endl;
        }
    }
    else
    {
        WriterEmpty write;
        sampler.generateUniform(nbSample, -1, write, seed);
    }
}

