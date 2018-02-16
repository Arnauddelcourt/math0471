// example of hybrid MPI/OpenMP program
//     run with 2 processes and 6 threads per process (ubuntu)
//         export OMP_NUM_THREADS=6
//         [ubuntu - openmpi]
//         mpirun -np 2 -cpus-per-rank 6 --bind-to core:overload-allowed  bin/fdtd_mpi
//         [windows - microsoft mpi]
//         mpiexec -np 2 bin\fdtd_mpi

#include "vtl.h"
#include "vtlSPoints.h"

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <cassert>

#include <mpi.h>
#include <omp.h>
#include <stdlib.h>

using namespace vtl;

int main(int argc, char *argv[])
{
    // Global grid parameters

    SPoints grid;
    grid.o = Vec3d(10.0, 10.0, 10.0); // origin
    Vec3d L(50.0, 60.0, 80.0); // box dimensions

    grid.np1 = Vec3i(0, 0, 0);       // first index
    grid.np2 = Vec3i(100, 120, 160); // last index

    //grid.np = grid.np2 - grid.np1 + 1; // nb of points
    grid.dx = L / (grid.np() - 1);  // compute spacing

    int nstepT = 20; // nb of time steps
    int saveResults = true;

    // MPI init
    MPI_Init(&argc, &argv);
    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // display some OpenMP/MPI info
    if (myid == 0)
    {
        std::cout << myid << ": we have " << numprocs << " processes\n";
        int nthreads = omp_get_max_threads();
        std::cout << myid << ": we have " << nthreads << " threads\n";
		char *omp_num_t = getenv("OMP_NUM_THREADS");
		
		std::cout << "OMP_NUM_THREADS="	
					<< (omp_num_t? omp_num_t : "undefined") << '\n';
    }

    // get my grid indices

    SPoints mygrid = grid.split(numprocs, myid);
    std::cout << myid << ": " << mygrid << '\n';

    // build extents vectors

    std::vector<SPoints> sgrids; // list of subgrids (used by rank0 for pvti file)
    if (myid == 0)
    {
        sgrids.resize(numprocs);
        for (int i = 0; i < numprocs; ++i)
            sgrids[i] = grid.split(numprocs, i);
    }

    // creation of dummy fields (over my subdomain)
    int mynbp = mygrid.nbp();

    std::cout << myid << ": " << mynbp << " points created\n";

    std::vector<double> scalarX(mynbp);
    std::vector<double> scalarY(mynbp);
    std::vector<double> scalarZ(mynbp);
    std::vector<double> vectorSIN(mynbp * 3); // vector field

    mygrid.scalars["scalar X"] = &scalarX;
    mygrid.scalars["scalar Y"] = &scalarY;
    mygrid.scalars["scalar Z"] = &scalarZ;
    mygrid.vectors["vector SIN"] = &vectorSIN;

    //MPI_Barrier(MPI_COMM_WORLD);

    // NOTE: the origin is the same for all subdomains

    // time step loop
    for (int nstep = 0; nstep < nstepT; ++nstep)
    {
        double time = double(nstep) / nstepT;

        //int idx = 0;
        int npz1 = mygrid.np1[2];
        int npz2 = mygrid.np2[2];
        Vec3i np = mygrid.np();

#pragma omp parallel for
        for (int k = npz1; k <= npz2; ++k)
        {
            double z = k * mygrid.dx[2] + mygrid.o[2];
            int npy1 = mygrid.np1[1];
            int npy2 = mygrid.np2[1];            
            for (int j = npy1; j <= npy2; ++j)
            {
                double y = j * mygrid.dx[1] + mygrid.o[1];
                int npx1 = mygrid.np1[0];
                int npx2 = mygrid.np2[0];                 
                for (int i = npx1; i <= npx2; ++i)
                {
                    double x = i * mygrid.dx[0] + mygrid.o[0];

                    int idx = (k-npz1) * (np[1] * np[0]) + (j-npy1) * np[0] + (i-npx1);
                    //idx++;

                    scalarX[idx] = x;
                    scalarY[idx] = y;
                    scalarZ[idx] = z;
                    vectorSIN[idx * 3 + 0] = sin(2 * M_PI * (((x - (mygrid.o[0] + L[0] / 2.)) / L[0]) + time));
                    vectorSIN[idx * 3 + 1] = cos(2 * M_PI * (((y - (mygrid.o[1] + L[1] / 2.)) / L[1]) + time));
                    vectorSIN[idx * 3 + 2] = sin(2 * M_PI * (((z - (mygrid.o[2] + L[2] / 2.)) / L[2]) + time));
                }
            }
        }
        if(saveResults)
        {
            // save results of the mpi process to disk
            export_spoints_XML("fdtd", nstep, grid, mygrid, UNZIPPED);
            export_spoints_XML("fdtdz", nstep, grid, mygrid, ZIPPED);

            if (myid == 0) // save main pvti file by rank0
            {
                export_spoints_XMLP("fdtd", nstep, grid, mygrid, sgrids, UNZIPPED);
                export_spoints_XMLP("fdtdz", nstep, grid, mygrid, sgrids, ZIPPED);
            }
        }
    }

    MPI_Finalize();

    return 0;
}
