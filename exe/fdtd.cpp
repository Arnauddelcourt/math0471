#include "cube.h"
#include "vtklite.h"

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <cassert>

int main(int argc, char *argv[])
{
    // creation of a cube of particles

    double o[3] = { 10.0, 10.0, 10.0 };
    double L[3] = { 50.0, 60.0, 80.0 };
    int np[3] = { 50, 60, 80 };
    int nstepT = 1;

    // creation of dummy fields

    int nbp = np[0]*np[1]*np[2];

    double dx[3];
    for(int i=0; i<3; ++i)
        dx[i] = L[i]/(np[i]-1);

    std::cout << nbp << " points created\n";
    std::vector<double> Ez(nbp);
    std::vector<double> Hy(nbp);
    std::vector<double> poynting(nbp*3);

    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;
    scalars["Ez"] = &Ez;
    scalars["Hy"]  = &Hy;
    vectors["poynting"] = &poynting;

    // time step loop
    for(int nstep=0; nstep<nstepT; ++nstep)
    {   
        int idx = 0;
        for(int k=0; k<np[2]; ++k)
        {
            double z = k*dx[2]+o[2]; 
            for(int j=0; j<np[1]; ++j)
            {
                double y = j*dx[1]+o[1];
                for(int i=0; i<np[0]; ++i)
                {
                    double x = i*dx[0]+o[0];
                    
                    Ez[idx] = x;
                    Hy[idx] = y;
                    poynting[idx*3+0] = x;
                    poynting[idx*3+1] = y;
                    poynting[idx*3+2] = z;
                    idx++;
                }
        
            }
        }

        // save results to disk
        paraview2("results_txt", nstep, o, dx, np, scalars, vectors, LEGACY_TXT);
        paraview2("results_bin", nstep, o, dx, np, scalars, vectors, LEGACY_BIN);
        paraview2("results_bin", nstep, o, dx, np, scalars, vectors, XML_BIN);
        paraview2("results_binz", nstep, o, dx, np, scalars, vectors, XML_BINZ);
    }

    return 0;
}
