#ifndef VTLSPOINTS_H
#define VTLSPOINTS_H

#include "vtl.h"
#include "vtlVec3.h"
#include <map>
#include <string>
#include <vector>


namespace vtl
{

// a structured grid

class VTL_API SPoints
{
  public:
    int id;     ///< rank of the grid
    Vec3d o;    ///< origin
    //Vec3d L;    ///< length of the domain along x,y,z
    Vec3i np1;  ///< starting indices
    Vec3i np2;  ///< ending indices
    //Vec3i np;   ///< nb of points along each direction
    Vec3d dx;   ///< spacing

    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;


    SPoints();
    SPoints split(int numprocs, int myid);
    Vec3d L() const { return Vec3d(np2-np1)*dx; }
    Vec3i np() const { return np2-np1+1; }
    int nbp() const { Vec3i a = np(); return a[0]*a[1]*a[2]; }
    friend std::ostream &operator<<(std::ostream &out, SPoints const &obj);   
};
}

#endif //VTLSPOINTS_H