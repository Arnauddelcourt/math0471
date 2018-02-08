FILE(GLOB SRCS *.h *.cpp)
ADD_EXECUTABLE(fdtd_mpi ${SRCS})

# MPI --
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(fdtd_mpi ${MPI_LIBRARIES})
# --

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(fdtd_mpi vtklite)
