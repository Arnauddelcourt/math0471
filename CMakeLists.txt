PROJECT(SPH CXX)
cmake_minimum_required(VERSION 2.7)

ADD_SUBDIRECTORY( src )
ADD_SUBDIRECTORY( exe/sph )
ADD_SUBDIRECTORY( exe/fdtd )
ADD_SUBDIRECTORY( exe/fdtd_mpi )

