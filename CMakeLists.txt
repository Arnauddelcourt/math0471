PROJECT(SPH CXX)
cmake_minimum_required(VERSION 2.7)

FILE(GLOB SRCS *.h *.cpp)
ADD_EXECUTABLE(sph ${SRCS})




