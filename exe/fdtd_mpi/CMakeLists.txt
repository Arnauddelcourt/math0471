FILE(GLOB SRCS *.h *.cpp)
ADD_EXECUTABLE(fdtd_mpi ${SRCS})

# MPI --
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(fdtd_mpi ${MPI_LIBRARIES})
IF( MINGW )
    add_definitions(-DMSMPI_NO_SAL)
ENDIF()


# OpenMP --
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    #MESSAGE ("OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()



INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(fdtd_mpi vtklite)
