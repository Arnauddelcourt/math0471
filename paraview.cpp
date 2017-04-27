#include "sph.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdint.h>
#include "swapbytes.h"


const int __one__ = 1;
const bool isCpuLittleEndian = 1 == *(char*)(&__one__); // CPU endianness

// this routine writes a double vector to file f.
//  the vector is converted to float (32bits) and to big endian format (required by the legacy VTK format)

void write_vector(std::ofstream &f, std::vector<double> const &pos, int nbp, int nj, bool ascii)
{
    assert(pos.size()==nbp*nj);
    if(ascii) 
    {
        for(int i=0; i<nbp; ++i)
        {
            for(int j=0;j<nj;++j)
                f << pos[nj*i+j] << " ";
            f << '\n';
        }
    }
    else
    {
        if(isCpuLittleEndian)
            for(int i=0; i<nbp; ++i)
            {
                // float+little endian => double should be converted to float, then swapped
                for(int j=0; j<nj; ++j)
                {
                    float fx = (float)pos[nj*i+j]; 
                    uint32_t x = swap_uint32(*(uint32_t*)&fx); 
                    f.write((char*)&x, sizeof(uint32_t));
                }
            }
        else
        {
            // double+bigendian => vector can be written as in memory
            //f.write(reinterpret_cast<char const*>(&(pos[0])), pos.size()*sizeof(double));
            // float+bigendian => vector should be converted to float
            for(int i=0; i<nbp; ++i)
                for(int j=0; j<nj; ++j)
                {
                    float fx = (float)pos[nj*i+j]; 
                    f.write((char*)&fx, sizeof(uint32_t));
                }
        }
    }
}

// export results to paraview (VTK polydata - legacy file fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]

void paraview(std::string const &filename, 
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors, bool ascii)
{
    //std::cout << "system is " << (isCpuLittleEndian? "little" : "big") << " endian\n";
    //bool ascii=false;

    int nbp = pos.size()/3;
    assert(pos.size()==nbp*3); // should be multiple of 3
    
    // build file name + stepno + vtk extension
    std::stringstream s; s << filename << std::setw(8) << std::setfill('0') << step << ".vtk";

    // open file
    std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    f << std::scientific;
    // header
    f << "# vtk DataFile Version 3.0\n";
    f << "file written by sph.exe\n";
    f << (ascii ? "ASCII\n" : "BINARY\n");
    f << "DATASET POLYDATA\n";

    // points
    f << "POINTS " << nbp << " float\n";
    write_vector(f, pos, nbp, 3, ascii);
    
    // vertices
    f << "VERTICES " << nbp << " " << 2*nbp << "\n";
    if(ascii) 
    {    
        for(int i=0; i<nbp; ++i)
            f << "1 " << i << '\n';
        f << '\n'; // empty line (required)
    }
    else
    {
        int32_t type = isCpuLittleEndian? swap_int32(1) : 1;
        for(int i=0; i<nbp; ++i)
        {
            int32_t ii = isCpuLittleEndian? swap_int32(i) : i;
            f.write((char*)&type, sizeof(int));
            f.write((char*)&ii, sizeof(int));
        }
    }

    // fields
    f << "POINT_DATA " << nbp << '\n';
    f << "FIELD FieldData " << scalars.size()+vectors.size() << '\n';

    // scalar fields
    std::map<std::string, std::vector<double> *>::const_iterator it=scalars.begin();
    for(; it!=scalars.end(); ++it)
    {
        assert(it->second->size()==nbp);
        f << it->first << " 1 " << nbp << " float\n";
        write_vector(f, *it->second, nbp, 1, ascii);
    }

    // vector fields
    it = vectors.begin();
    for(; it!=vectors.end(); ++it)
    {
        assert(it->second->size()==3*nbp);
        f << it->first << " 3 " << nbp << " float\n";
        write_vector(f, *it->second, nbp, 3, ascii);
    }
    f.close();
}
