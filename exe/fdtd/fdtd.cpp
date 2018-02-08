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

    double o[3] = {10.0, 10.0, 10.0};
    double L[3] = {50.0, 60.0, 80.0};
    int np[3] = {50, 60, 80};
    int nstepT = 20;

    // creation of dummy fields
    int nbp = np[0] * np[1] * np[2];

    double dx[3];
    for (int i = 0; i < 3; ++i)
        dx[i] = L[i] / (np[i] - 1);

    std::cout << nbp << " points created\n";
    std::vector<double> Ez(nbp);
    std::vector<double> Hy(nbp);
    std::vector<double> poynting(nbp * 3); // vector field

    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;
    scalars["Ez"] = &Ez;
    scalars["Hy"] = &Hy;
    vectors["poynting"] = &poynting;

    // time step loop
    for (int nstep = 0; nstep < nstepT; ++nstep)
    {
        double time = double(nstep) / nstepT;

        int idx = 0;
        for (int k = 0; k < np[2]; ++k)
        {
            double z = k * dx[2] + o[2];
            for (int j = 0; j < np[1]; ++j)
            {
                double y = j * dx[1] + o[1];
                for (int i = 0; i < np[0]; ++i)
                {
                    double x = i * dx[0] + o[0];

                    Ez[idx] = x;
                    Hy[idx] = y;
                    poynting[idx * 3 + 0] = sin(2*M_PI*(((x-(o[0]+L[0]/2.))/L[0])+time)) ;
                    poynting[idx * 3 + 1] = cos(2*M_PI*(((y-(o[1]+L[1]/2.))/L[1])+time)) ;
                    poynting[idx * 3 + 2] = sin(2*M_PI*(((z-(o[2]+L[2]/2.))/L[2])+time)) ;
                    idx++;
                }
            }
        }

        // save results to disk
        export_spoints("fdtd_legacy_ascii", nstep, o, dx, np, scalars, vectors, LEGACY_TXT);
        export_spoints("fdtd_legacy_bin", nstep, o, dx, np, scalars, vectors, LEGACY_BIN);
        export_spoints("fdtd_xml_apprawbin", nstep, o, dx, np, scalars, vectors, XML_BIN);
        export_spoints("fdtd_xml_apprawbinz", nstep, o, dx, np, scalars, vectors, XML_BINZ);
    }

    return 0;
}
