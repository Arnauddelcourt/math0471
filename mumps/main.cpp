// solves a Laplacian over a cube with mumps

#include "vtl.h"
#include "vtlSPoints.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>

using namespace vtl;

// from c_example.c ------
#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
// -----------------------

int main(int argc, char *argv[])
{
    bool matlab = false; // save matrix to file [debug]

    int error = 0;
    MUMPS_INT ierr = MPI_Init(&argc, &argv);
    MUMPS_INT myid;

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // -----------------------------------------------
    // init MUMPS
    // -----------------------------------------------

    DMUMPS_STRUC_C id;
    id.comm_fortran = USE_COMM_WORLD;
    id.par = 1; // 1=host involved in factorization phase
    id.sym = 0; // 0=unsymmetric
    id.job = JOB_INIT;
    std::cout << "Init MUMPS package.\n";
    dmumps_c(&id);
    if (id.infog[0] < 0)
    {
        printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
               myid, id.infog[0], id.infog[1]);
        error = 1;
    }
    else
    {
        std::cout << "OK!\n";
    }
    // -------------------------------------------------

    SPoints grid;



        // setup grid

        grid.o = Vec3d(10.0, 10.0, 10.0); // origin
        Vec3d L(20.0, 30.0, 40.0);        // box dimensions

        grid.np1 = Vec3i(0, 0, 0);    // first index
        grid.np2 = Vec3i(20, 30, 40); // last index

        grid.dx = L / (grid.np() - 1); // compute spacing

        // creation of dummy fields
        int nbp = grid.nbp();

     
        std::cout << nbp << " points created\n";
        std::cout << grid;

        std::vector<MUMPS_INT> irn;
        std::vector<MUMPS_INT> jcn;
        std::vector<double> A;
        std::vector<double> rhs; //(grid.nbp());
        grid.scalars["Temp"] = &rhs;

        // build matrix & rhs

        if (myid == 0)
        {         
            
        int npz1 = grid.np1[2];
        int npz2 = grid.np2[2];
        int npy1 = grid.np1[1];
        int npy2 = grid.np2[1];
        int npx1 = grid.np1[0];
        int npx2 = grid.np2[0];
        Vec3i np = grid.np();

        auto loc = [=](int i, int j, int k) { return (k - npz1) * (np[1] * np[0]) + (j - npy1) * np[0] + (i - npx1) + 1; }; // +1 (fortran)

        double dx2 = grid.dx[0] * grid.dx[0];
        double dy2 = grid.dx[1] * grid.dx[1];
        double dz2 = grid.dx[2] * grid.dx[2];

        for (int k = npz1; k <= npz2; ++k)
        {
            double z = k * grid.dx[2] + grid.o[2];
            for (int j = npy1; j <= npy2; ++j)
            {
                double y = j * grid.dx[1] + grid.o[1];
                for (int i = npx1; i <= npx2; ++i)
                {
                    double x = i * grid.dx[0] + grid.o[0];

                    int id = loc(i, j, k);

                    // BCs ===========================================

                    if (k == npz1) // impose dirichlet first
                    {
                        // dirichlet
                        irn.push_back(id);
                        jcn.push_back(id);
                        A.push_back(1.0);
                        //irn.push_back(id);
                        //jcn.push_back(loc(i, j, k + 1));
                        //A.push_back(-1.0);
                    }
                    else if (k == npz2)
                    {
                        // dirichlet
                        irn.push_back(id);
                        jcn.push_back(id);
                        A.push_back(1.0);
                        //irn.push_back(id);
                        //jcn.push_back(loc(i, j, k - 1));
                        //A.push_back(-1.0);
                    }
                    else if (j == npy1)
                    {
                        // neumann
                        irn.push_back(id);
                        jcn.push_back(id);
                        A.push_back(1.0);
                        irn.push_back(id);
                        jcn.push_back(loc(i, j + 1, k));
                        A.push_back(-1.0);
                    }
                    else if (j == npy2)
                    {
                        // neumann
                        irn.push_back(id);
                        jcn.push_back(id);
                        A.push_back(1.0);
                        irn.push_back(id);
                        jcn.push_back(loc(i, j - 1, k));
                        A.push_back(-1.0);
                    }
                    else if (i == npx1)
                    {
                        // neumann
                        irn.push_back(id);
                        jcn.push_back(id);
                        A.push_back(1.0);
                        irn.push_back(id);
                        jcn.push_back(loc(i + 1, j, k));
                        A.push_back(-1.0);
                    }
                    else if (i == npx2)
                    {
                        // neumann
                        irn.push_back(id);
                        jcn.push_back(loc(i - 1, j, k));
                        A.push_back(1.0);
                        irn.push_back(id);
                        jcn.push_back(loc(i + 1, j, k));
                        A.push_back(-1.0);
                    }
                    else
                    {

                        if ((i != npx1) && (i != npx2) && (j != npy1) && (j != npy2) && (k != npz1) && (k != npz2))
                        {
                            irn.push_back(id);
                            jcn.push_back(id);
                            A.push_back(-2.0 * (1.0 / dx2 + 1.0 / dy2 + 1.0 / dz2));
                        }
                        // x
                        if (i != npx1)
                        {
                            irn.push_back(id);
                            jcn.push_back(loc(i - 1, j, k));
                            A.push_back(1.0 / dx2);
                        }
                        if (i != npx2)
                        {
                            irn.push_back(id);
                            jcn.push_back(loc(i + 1, j, k));
                            A.push_back(1.0 / dx2);
                        }
                        // y
                        if (j != npy1)
                        {
                            irn.push_back(id);
                            jcn.push_back(loc(i, j - 1, k));
                            A.push_back(1.0 / dy2);
                        }
                        if (j != npy2)
                        {
                            irn.push_back(id);
                            jcn.push_back(loc(i, j + 1, k));
                            A.push_back(1.0 / dy2);
                        }

                        // y
                        if (k != npz1)
                        {
                            irn.push_back(id);
                            jcn.push_back(loc(i, j, k - 1));
                            A.push_back(1.0 / dz2);
                        }
                        if (k != npz2)
                        {
                            irn.push_back(id);
                            jcn.push_back(loc(i, j, k + 1));
                            A.push_back(1.0 / dz2);
                        }
                    }

                    // rhs
                    if (k == npz1)
                    {
                        rhs.push_back(1.0);
                    }
                    else if (k == npz2)
                    {
                        rhs.push_back(2.0);
                    }
                    else
                    {
                        rhs.push_back(0.0);
                    }
                    //std::cout << "rhs.size()=" << rhs.size() << '\n';
                }
            }
        }
/*
        std::ofstream f("matrix.m");
        if (matlab)
        {
            std::cout << "saving matrix to file...\n";

            for (size_t i = 0; i < A.size(); ++i)
            {
                f << "A(" << irn[i] << "," << jcn[i] << ")=" << A[i] << ";\n";
            }
            for (size_t i = 0; i < rhs.size(); ++i)
            {
                f << "rhs(" << i + 1 << ")=" << rhs[i] << ";\n";
            }
        }
*/
        std::cout << "nnz=" << A.size() << '\n';

        /* Define the problem on the host */
        id.n = rhs.size();
        id.nnz = A.size();
        id.irn = &irn[0];
        id.jcn = &jcn[0];
        id.a = &A[0];
        id.rhs = &rhs[0];
    }

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
    /* No outputs */
    id.ICNTL(1) = -1; // stream for error messages [def=6]
    id.ICNTL(2) = -1; // stream for diag printing, statistics, warnings [def=0]
    id.ICNTL(3) = -1; // stream for global information [def=6]
    id.ICNTL(4) = 0;  // level of printing [def=2]
    // id.ICNTL(5)   // matrix input format
    // id.ICNTL(6)   // permutation/scaling
    // id.ICNTL(7)   // ordering
    // id.ICNTL(8)   // scaling strategy [def=auto]
    // id.ICNTL(9)   // use A or A^T [def=A]
    // id.ICNTL(10)  // iterative refinement [def=0=disabled]
    // id.ICNTL(11)  // compute statistics on error [def=0=disabled]
    // id.ICNTL(12)  // ordering strategy for sym matrices [def=0=auto]
    // id.ICNTL(13)  // parallelism of root node (scalapack) [def=0=parallel with scalapack]
    // id.ICNTL(14)  // % incr of working space [def=20=20%]
    // id.ICNTL(15-17)  // NOT USED
    // id.ICNTL(18)  // distributed input matrix [def=0=centralized]
    // id.ICNTL(19)  // Schur complement [def=0=no schur cplt]
    // id.ICNTL(20)  // format of rhs [def=0=dense]
    // id.ICNTL(21)  // distribution of solution vectors [def=0=centralized]
    // id.ICNTL(22)  // out-of-core [def=0=in-core]
    // id.ICNTL(23)  // max memory [def=0=estimate]
    // id.ICNTL(24)  // null pivot detectio [def=0=disabled]
    // id.ICNTL(25)  // solution for deficiant matrix [def=0=1 sol is returned]
    // id.ICNTL(26)  // [see schur cplt]
    // id.ICNTL(27)  // blocking size for multiple rhs [def=-32=auto]
    // id.ICNTL(28)  // parallel ordering [def=0=auto]
    // id.ICNTL(29)  // parallel ordering method (if scotch/parmetis) [def=0=auto]
    // id.ICNTL(30)  // compute some A^-1 entries
    // id.ICNTL(31)  // keep factors [def=0=yes]
    // id.ICNTL(32)  // forward elimination during factorization [def=0=disabled]
    // id.ICNTL(33)  // compute det(A)
    // id.ICNTL(34)  // NOT USED
    // id.ICNTL(35)  // BLR factorization (def=0=disabled)

    /* Call the MUMPS package (analyse, factorization and solve). */
    std::cout << "Call the MUMPS package (analyse, factorization and solve).\n";
    id.job = 6;
    dmumps_c(&id);
    if (id.infog[0] < 0)
    {
        printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
               myid, id.infog[0], id.infog[1]);
        error = 1;
    }
    else
    {
        std::cout << "OK!\n";
    }

    /* Terminate instance. */
    std::cout << "Terminate instance.\n";
    id.job = JOB_END;
    dmumps_c(&id);
    if (myid == 0)
    {
        if (!error)
        {
            printf("Solution is : (%8.2f  %8.2f)\n", id.rhs[0], id.rhs[1]);
        }
        else
        {
            printf("An error has occured, please check error code returned by MUMPS.\n");
        }
    }

    if (myid == 0)
    {
        if (matlab)
        {
            /*
            std::cout << "saving solution to file...\n";
            for (size_t i = 0; i < id.n; ++i)
            {
                f << "sol(" << i + 1 << ")=" << id.rhs[i] << ";\n";
            }
            f << "[A\\rhs' sol']\n";
            f << "error = norm(A\\rhs'-sol')\n";
            f.close();
            */
        }

        // save results to disk
        export_spoints_XML("laplace", 0, grid, grid, Zip::ZIPPED);
    }

    ierr = MPI_Finalize();

    return 0;
}
