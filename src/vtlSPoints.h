#ifndef VTLSPOINTS_H
#define VTLSPOINTS_H

#include "vtl.h"
#include "vtlVec3.h"

namespace vtl
{

// a structured grid

class VTL_API SPoints
{
  public:
    Vec3d o;
    Vec3d L;
    Vec3i np;
    Vec3d dx;

    SPoints();
};
}

#endif //VTLSPOINTS_H