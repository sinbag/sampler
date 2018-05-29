#include "write-eps.h"

#include <fstream>

// export sites to an EPS image
template <typename T>
void write_eps(std::string filename, const std::vector<T> &points, int nDims, double radius, double scale) {

    if(filename.compare(filename.size()-4, 4,".eps") != 0){
        filename.erase(filename.end()-4, filename.end());
        filename += ".eps";
    }


    std::ofstream os;
    os.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);

    radius /= scale;

    os << "%!PS-Adobe-3.1 EPSF-3.0\n";
    os << "%%HiResBoundingBox: " << -radius << " " << -radius << " " << scale+radius << " " << scale+radius << "\n";
    os << "%%BoundingBox: " << -radius << " " << -radius << " " << scale+radius << " " << scale+radius << "\n";
    os << "%%CropBox: " << -radius << " " << -radius << " " << scale+radius << " " << scale+radius << "\n";
    os << "/radius { " << radius << " } def\n";
    os << "/p { radius 0 360 arc closepath fill stroke } def\n";
    os << "gsave " << scale << " " << scale << " scale\n";

    os << "0 0 0 setrgbcolor\n";

    if(nDims == 1){
        int npoints = points.size();
        for (unsigned int i = 0; i < npoints; ++i) {
            os << points[i]+0.5 << " " << 0.5 << " p\n";
        }
    }
    else{
        int npoints = points.size() * 0.5;
        for (unsigned int i = 0; i < npoints; ++i) {
            os << points[2*i] << " " << points[2*i+1] << " p\n";
        }
    }
    os << "grestore\n";
    os.close();

    //  return true;
}

template void write_eps(std::string filename, const std::vector<double> &points, int nDims, double radius, double scale);
template void write_eps(std::string filename, const std::vector<float> &points, int nDims, double radius, double scale);

// export sites to an EPS image
void write_eps_colored(std::string filename, const std::vector<double> &points, int nDims, double radius, double scale) {

    if(filename.compare(filename.size()-4, 4,".eps") != 0){
        filename.erase(filename.end()-4, filename.end());
        filename += ".eps";
    }


    std::ofstream os;
    os.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);

    radius /= scale;

    os << "%!PS-Adobe-3.1 EPSF-3.0\n";
    os << "%%HiResBoundingBox: " << -radius << " " << -radius << " " << scale+radius << " " << scale+radius << "\n";
    os << "%%BoundingBox: " << -radius << " " << -radius << " " << scale+radius << " " << scale+radius << "\n";
    os << "%%CropBox: " << -radius << " " << -radius << " " << scale+radius << " " << scale+radius << "\n";
    os << "/radius { " << radius << " } def\n";
    os << "/p { radius 0 360 arc closepath fill stroke } def\n";
    os << "gsave " << scale << " " << scale << " scale\n";

    os << "0 0 0 setrgbcolor\n";

    if(nDims == 1){
        int npoints = points.size();
        for (unsigned int i = 0; i < npoints; ++i) {
            os << points[i]+0.5 << " " << 0.5 << " p\n";
        }
    }
    else{
        int npoints = points.size() * 0.5;
        for (unsigned int i = 0; i < npoints; ++i) {
            os << points[2*i] << " " << points[2*i+1] << " p\n";
        }
    }
    os << "grestore\n";
    os <<"1 0 0 setrgbcolor\n";
    os.close();

    //  return true;
}

// export sites to an TXT image
template <typename T>
void write_txt(std::string filename, const std::vector<T> &points, int nDims) {
    
    if(filename.compare(filename.size()-4, 4,".txt") != 0){
        filename.erase(filename.end()-4, filename.end());
        filename += ".txt";
    }
    
    std::ofstream os;
    os.open(filename.c_str());
    
    int npoints = points.size() * 0.5;
    for (unsigned int i = 0; i < npoints; ++i) {
        os << points[2*i] << " " << points[2*i+1] << " \n";
    }
    os.close();
}

template void write_txt(std::string filename, const std::vector<double> &points, int nDims);
template void write_txt(std::string filename, const std::vector<float> &points, int nDims);


