#include "sph.h"
#include <fstream>

void paraview(std::string const &filename,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors)
{
    int nbp = pos.size()/3;

    std::ofstream f(filename+".vtk");
    f << "# vtk DataFile Version 3.0\n";
    f << "file written by sph.exe\n";
    f << "ASCII\n";
    f << "DATASET POLYDATA\n";
    f << "POINTS " << nbp << " float\n";
    for(int i=0; i<nbp; ++i)
        f << pos[3*i+0] << " " << pos[3*i+1] << " " << pos[3*i+2] << '\n';
    f << "VERTICES " << nbp << " " << 2*nbp << "\n";
    for(int i=0; i<nbp; ++i)
        f << "1 " << i << '\n';
    f << '\n';
    f << "POINT_DATA " << nbp << '\n';
    f << "FIELD FieldData " << scalars.size()+vectors.size() << '\n';
    std::map<std::string, std::vector<double> *>::const_iterator it=scalars.begin();
    for(; it!=scalars.end(); ++it)
    {
        f << it->first << " 1 " << nbp << " float\n";
        for(int i=0; i<nbp; ++i)
            f << (*it->second)[i] << '\n';
    }
    it = vectors.begin();
    for(; it!=vectors.end(); ++it)
    {
        f << it->first << " 3 " << nbp << " float\n";
        for(int i=0; i<nbp; ++i)
            f << (*it->second)[3*i+0] << " " << (*it->second)[3*i+1] << " " << (*it->second)[3*i+2] << '\n';
    }
    f.close();
}