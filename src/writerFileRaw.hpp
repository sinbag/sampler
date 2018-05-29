#ifndef WRITER_FILE_RAW_H
#define WRITER_FILE_RAW_H


#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "tileState.hpp"
#include "vector.hpp"

class WriterFileRaw
{
    private:
    std::vector<float> m_pointset;
        std::string m_filename;
        std::ofstream m_file;
        float m_spaceSize;

    public:
        WriterFileRaw(const std::string& fn);
        ~WriterFileRaw();

        void open(const float& spaceSize=1.);

    public:
        inline void sample(const Point& sample, const TileState& tilestate)
        {
            //m_file << sample.x()/m_spaceSize+0.5 << "\t" << sample.y()/m_spaceSize+0.5 << std::endl;
            m_pointset.push_back(sample.x()/m_spaceSize+0.5);
            m_pointset.push_back(sample.y()/m_spaceSize+0.5);
        }

        void clear(const float& spaceSize=1.);
        const std::vector<float>& pts() const { return m_pointset; }
};

#endif

