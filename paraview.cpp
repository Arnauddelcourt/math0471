#include "sph.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
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
//   ascii:   'true' for ASCII format, 'false' for binary

void paraview(std::string const &filename, 
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors, 
              bool ascii)
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

size_t write_vectorXML(std::ofstream &f, std::vector<double> const &pos, int nbp, int nj)
{
    assert(pos.size()==nbp*nj);
    size_t written=0;
    // data block size
    uint32_t sz = nbp*nj*sizeof(uint32_t);
    f.write((char*)&sz, sizeof(uint32_t));
    written+=sizeof(uint32_t);
    // data
    for(int i=0; i<nbp; ++i)
        for(int j=0; j<nj; ++j)
        {
            float fx = (float)pos[nj*i+j]; 
            f.write((char*)&fx, sizeof(uint32_t));
        }
    written+=sz;
    return written;
}

// export results to paraview (VTK polydata - XML fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]

void paraviewXML(std::string const &filename, 
                 int step,
                 std::vector<double> const &pos,
                 std::map<std::string, std::vector<double> *> const &scalars,
                 std::map<std::string, std::vector<double> *> const &vectors, 
                 bool ascii)
{
    //std::cout << "system is " << (isCpuLittleEndian? "little" : "big") << " endian\n";
    //bool ascii=false;

    int nbp = pos.size()/3;
    assert(pos.size()==nbp*3); // should be multiple of 3
    
    // build file name + stepno + vtk extension
    std::stringstream s; s << filename << std::setw(8) << std::setfill('0') << step << ".vtp";
    std::stringstream s2; s2 << filename << std::setw(8) << std::setfill('0') << step << ".vtp.tmp";

    // open file
    std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    std::ofstream f2(s2.str().c_str(), std::ios::binary | std::ios::out); // temp binary file
    f << std::scientific;

    size_t offset = 0;
    // header
    f << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <PolyData>\n";
    f << "    <Piece NumberOfPoints=\"" << nbp << "\" ";
    f << "NumberOfVerts=\"" << nbp << "\" ";
    f << "NumberOfLines=\"0\" ";
    f << "NumberOfStrips=\"0\" ";
    f << "NumberOfPolys=\"0\">\n";

    // ------------------------------------------------------------------------------------
    f << "      <PointData>\n";
    // scalar fields
    std::map<std::string, std::vector<double> *>::const_iterator it=scalars.begin();
    for(; it!=scalars.end(); ++it)
    {
        assert(it->second->size()==nbp);
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << it->first << "\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *it->second, nbp, 1);
    }
    // vector fields
    it = vectors.begin();
    for(; it!=vectors.end(); ++it)
    {
        assert(it->second->size()==3*nbp);
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << it->first << "\" ";
        f << " NumberOfComponents=\"3\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *it->second, nbp, 3);
    }
    f << "      </PointData>\n";

    // ------------------------------------------------------------------------------------
    f << "      <CellData>\n";
    f << "      </CellData>\n";  

    // ------------------------------------------------------------------------------------
    f << "      <Points>\n";
    f << "        <DataArray type=\"Float32\" ";
    f << " Name=\"Points\" ";
    f << " NumberOfComponents=\"3\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";   
    offset += write_vectorXML(f2, pos, nbp, 3);
    f << "      </Points>\n";
    // ------------------------------------------------------------------------------------
    f << "      <Verts>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"" << nbp-1 << "\" ";
    f << " offset=\"" << offset << "\" />\n";
    // data block size
    uint32_t sz = nbp*sizeof(int);
    f2.write((char*)&sz, sizeof(uint32_t));
    offset+=sizeof(uint32_t);
    // data
    for(int i=0; i<nbp; ++i)
        f2.write((char*)&i, sizeof(int));
    offset += sz;

    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"1\" ";
    f << " RangeMax=\"" << nbp << "\" ";
    f << " offset=\"" << offset << "\" />\n";
    // data block size
    sz = nbp*sizeof(int);
    f2.write((char*)&sz, sizeof(uint32_t));
    offset+=sizeof(uint32_t);
    // data
    for(int i=1; i<nbp+1; ++i)
        f2.write((char*)&i, sizeof(int));
    offset += sz;

    f << "      </Verts>\n";


    // ------------------------------------------------------------------------------------
    f << "      <Lines>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n"; 
    sz = 0;
    f2.write((char*)&sz, sizeof(uint32_t));
    offset+=sizeof(uint32_t);      
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n"; 
    sz = 0;
    f2.write((char*)&sz, sizeof(uint32_t));
    offset+=sizeof(uint32_t);
    f << "      </Lines>\n";

    // ------------------------------------------------------------------------------------
    f << "      <Strips>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";   
    sz = 0;
    f2.write((char*)&sz, sizeof(uint32_t));
    offset+=sizeof(uint32_t);
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n"; 
    sz = 0;
    f2.write((char*)&sz, sizeof(uint32_t));
    offset+=sizeof(uint32_t);
    f << "      </Strips>\n";

    // ------------------------------------------------------------------------------------
    f << "      <Polys>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";  
    sz = 0;
    f2.write((char*)&sz, sizeof(uint32_t));
    offset+=sizeof(uint32_t); 
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    sz = 0;
    f2.write((char*)&sz, sizeof(uint32_t));
    offset+=sizeof(uint32_t); 
    f << "      </Polys>\n";


    f2.close();

    // ------------------------------------------------------------------------------------
    f << "    </Piece>\n";
    f << "  </PolyData>\n";
    // ------------------------------------------------------------------------------------
    f << "  <AppendedData encoding=\"raw\">\n";
    f << "    _";

    std::ifstream f3(s2.str().c_str(), std::ios::binary | std::ios::in);
    f << f3.rdbuf();
    f3.close();

    std::remove(s2.str().c_str());

    // binary here
    f << "  </AppendedData>\n";
    f << "</VTKFile>\n";

    f.close();
}
