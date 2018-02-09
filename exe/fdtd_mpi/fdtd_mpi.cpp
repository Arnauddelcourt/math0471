#include "vtklite.h"

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <cassert>

#include <mpi.h>

void splitgrid(int np[3], int numprocs, int myid, int myn1[3], int myn2[3])
{
    int sdir = 0; // split direction

    // calculate nb of points in each subdomain => mynp[]
    int mynp;
    int nfloor = np[sdir] / numprocs;
    int rem = np[sdir] % numprocs;
    assert(np[sdir] == nfloor * numprocs + rem);

    // build the point idx separating each subdomain
    std::vector<int> nps(numprocs + 1);
    nps[0] = 0;
    for (int i = 0; i < numprocs+1; ++i)
    {
        nps[i + 1] = nps[i] + nfloor;
        if (i < rem)
            nps[i + 1]++;
    }
    /*
    if(myid==0) // DEBUG
    {
        std::cout << "nps=";
        for (int i = 0; i < numprocs+1; ++i)
            std::cout << nps[i] << " ";
        std::cout << "\n";
    }
    */
    for (int i = 0; i < 3; ++i)
    {
        myn1[i] = 0;
        myn2[i] = np[i]-1;
    }
    myn1[sdir] = nps[myid];
    myn2[sdir] = nps[myid + 1]-1; 

    // add extra layers from direct neighborhood
    if(myid!=0)
        myn1[sdir]-=1;
}

int main(int argc, char *argv[])
{
    // Global grid parameters
    double o[3] = {10.0, 10.0, 10.0};
    double L[3] = {50.0, 60.0, 80.0};
    int np[3] = {51, 61, 81}; // nb of points in each direction - idx starting from 0 to np[i]-1
    int nstepT = 20;

    // compute spacing
    double dx[3];
    for (int i = 0; i < 3; ++i)
        dx[i] = L[i] / (np[i] - 1);

    // MPI init
    MPI_Init(&argc, &argv);
    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0)
        std::cout << myid << ": we have " << numprocs << " processes\n";

    // get my grid indices
    int myn1[3], myn2[3];
    splitgrid(np, numprocs, myid, myn1, myn2);
    std::cout << myid << ": [" 
    << myn1[0] << '-' << myn2[0] << ", " 
    << myn1[1] << '-' << myn2[1] << ", " 
    << myn1[2] << '-' << myn2[2] << "]";
    
    int mynp[3];
    for (int i = 0; i < 3; ++i)
        mynp[i] = myn2[i] - myn1[i] +1;

    std::cout << " = [" << mynp[0] << ", " << mynp[1] << ", " << mynp[2] << "]\n"; 
    // build extents vectors
    std::vector<std::vector<int> > extents(numprocs);
    if (myid == 0)
    {
        for (int i = 0; i < numprocs; ++i)
        {
            int myn1[3], myn2[3];
            splitgrid(np, numprocs, i, myn1, myn2);
            extents[i].resize(6);
            extents[i][0] = myn1[0];
            extents[i][1] = myn2[0];
            extents[i][2] = myn1[1];
            extents[i][3] = myn2[1];
            extents[i][4] = myn1[2];
            extents[i][5] = myn2[2];
        }
    }

    // creation of dummy fields (over my subdomain)
    int mynbp = mynp[0] * mynp[1] * mynp[2];

    std::cout << myid << ": " << mynbp << " points created\n";
    std::vector<double> Ez(mynbp);
    std::vector<double> Hy(mynbp);
    std::vector<double> poynting(mynbp * 3); // vector field

    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;
    scalars["Ez"] = &Ez;
    scalars["Hy"] = &Hy;
    vectors["poynting"] = &poynting;

    // NOTE: the origin is the same for all subdomains

    // time step loop
    for (int nstep = 0; nstep < nstepT; ++nstep)
    {
        double time = double(nstep) / nstepT;

        int idx = 0;
        for (int k = 0; k < mynp[2]; ++k)
        {
            double z = (myn1[2] + k) * dx[2] + o[2];
            for (int j = 0; j < mynp[1]; ++j)
            {
                double y = (myn1[1] + j) * dx[1] + o[1];
                for (int i = 0; i < mynp[0]; ++i)
                {
                    double x = (myn1[0] + i) * dx[0] + o[0];

                    Ez[idx] = x;
                    Hy[idx] = y;
                    poynting[idx * 3 + 0] = sin(2 * M_PI * (((x - (o[0] + L[0] / 2.)) / L[0]) + time));
                    poynting[idx * 3 + 1] = cos(2 * M_PI * (((y - (o[1] + L[1] / 2.)) / L[1]) + time));
                    poynting[idx * 3 + 2] = sin(2 * M_PI * (((z - (o[2] + L[2] / 2.)) / L[2]) + time));
                    idx++;
                }
            }
        }

        // save results to disk
        //export_spoints("fdtd_legacy_ascii", nstep, o, dx, np, scalars, vectors, LEGACY_TXT, myid, myn1, myn2);
        //export_spoints("fdtd_legacy_bin", nstep, o, dx, np, scalars, vectors, LEGACY_BIN, myid, myn1, myn2);
        export_spoints("fdtd_xml_apprawbin", nstep, o, dx, np, scalars, vectors, XML_BIN, myid, myn1, myn2);
        export_spoints("fdtd_xml_apprawbinz", nstep, o, dx, np, scalars, vectors, XML_BINZ, myid, myn1, myn2);

        if (myid == 0)
        {
            export_spoints_XMLP("fdtd_xml_apprawbin", nstep, o, dx, np, scalars, vectors, true, false, extents);
            export_spoints_XMLP("fdtd_xml_apprawbinz", nstep, o, dx, np, scalars, vectors, true, true, extents);
        }
    }

    MPI_Finalize();
    return 0;
}
