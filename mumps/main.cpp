// solves a Laplacian over a cube with mumps

#include "vtl.h"
#include "vtlSPoints.h"
#include "laplace.h"

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
#define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation

// -----------------------

int get_my_rank()
{
    int myid;
    int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    return myid;
}

void check_MUMPS(DMUMPS_STRUC_C &id)
{
    if (id.infog[0] < 0)
    {
        std::cout << "[" << get_my_rank() << "] MUMPS Error:\n";
        std::cout << "\tINFOG(1)=" << id.infog[0] << '\n';
        std::cout << "\tINFOG(2)=" << id.infog[1] << std::endl;
    }
}

void init_MUMPS(DMUMPS_STRUC_C &id)
{
    id.comm_fortran = -987654; //USE_COMM_WORLD;
    id.par = 1;                // 1=host involved in factorization phase
    id.sym = 0;                // 0=unsymmetric
    id.job = -1;
    std::cout << "[" << get_my_rank() << "] Init MUMPS package." << std::endl;
    dmumps_c(&id);
    check_MUMPS(id);
}

void end_MUMPS(DMUMPS_STRUC_C &id)
{
    id.job = -2;
    std::cout << "[" << get_my_rank() << "] Terminate MUMPS instance." << std::endl;
    dmumps_c(&id);
    check_MUMPS(id);
}

void solve_MUMPS(DMUMPS_STRUC_C &id)
{
    /*
    id.ICNTL(1) = -1; // stream for error messages [def=6]
    id.ICNTL(2) = -1; // stream for diag printing, statistics, warnings [def=0]
    id.ICNTL(3) = -1; // stream for global information [def=6]
    id.ICNTL(4) = 0;  // level of printing [def=2]
    */
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

    std::cout << "[" << get_my_rank()
              << "] Call the MUMPS package (analyse, factorization and solve)." << std::endl;
    id.job = 6;
    dmumps_c(&id);
    check_MUMPS(id);
}

void host_work(DMUMPS_STRUC_C &id)
{
    bool matlab = false; // save matrix to file [debug]

    SPoints grid;

    // setup grid

    grid.o = Vec3d(10.0, 10.0, 10.0); // origin
    Vec3d L(20.0, 30.0, 40.0);        // box dimensions

    grid.np1 = Vec3i(0, 0, 0);    // first index
    grid.np2 = Vec3i(20, 20, 20); // last index

    grid.dx = L / (grid.np() - 1); // compute spacing

    // creation of dummy fields
    int nbp = grid.nbp();

    std::cout << nbp << " points created\n";
    std::cout << grid;

    std::vector<MUMPS_INT> irn;
    std::vector<MUMPS_INT> jcn;
    std::vector<double> A;
    std::vector<double> rhs;
    grid.scalars["Temp"] = &rhs;

    fill_system(grid, irn, jcn, A, rhs);
    if(matlab)
    {
        save_matrix("A", irn, jcn, A);
        save_vector("rhs", rhs);
    }
    std::cout << "nnz=" << A.size() << '\n';

    // Define the problem on the host
    id.n = rhs.size();
    id.nnz = A.size();
    id.irn = &irn[0];
    id.jcn = &jcn[0];
    id.a = &A[0];
    id.rhs = &rhs[0];

    solve_MUMPS(id);

    if (matlab)
        save_vector("sol", rhs);

    // save results to disk
    export_spoints_XML("laplace", 0, grid, grid, Zip::ZIPPED);
}

void slave_work(DMUMPS_STRUC_C &id)
{
    solve_MUMPS(id);  
}

int main(int argc, char *argv[])
{
    // initialise MUMPS/MPI
    MPI_Init(&argc, &argv);
    DMUMPS_STRUC_C id;
    init_MUMPS(id);

    // split work among processes
    if (get_my_rank() == 0)
        host_work(id);
    else
        slave_work(id);

    // finalise MUMPS/MPI
    end_MUMPS(id);
    MPI_Finalize();

    return 0;
}
