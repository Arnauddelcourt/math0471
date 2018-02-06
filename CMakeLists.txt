PROJECT(SPH CXX)
cmake_minimum_required(VERSION 2.7)

FIND_PACKAGE(ZLIB)
IF(ZLIB_FOUND)
    add_definitions(-DUSE_ZLIB)
    INCLUDE_DIRECTORIES( ${ZLIB_INCLUDE_DIR} )
ENDIF()

#FILE(GLOB SRCS *.h *.cpp)

SET(SRCS1 cube.h cube.cpp paraview.h paraview.cpp swapbytes.h)

ADD_EXECUTABLE(sph ${SRCS1} sph.h sph.cpp)
ADD_EXECUTABLE(fdtd ${SRCS1} fdtd.h fdtd.cpp)
IF(ZLIB_FOUND)
    TARGET_LINK_LIBRARIES(sph ${ZLIB_LIBRARY} )
    TARGET_LINK_LIBRARIES(fdtd ${ZLIB_LIBRARY} )
ENDIF()



