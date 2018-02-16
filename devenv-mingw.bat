@echo off

echo setting MinGW environment

set INCLUDE=c:\local\include
set INCLUDE=%INCLUDE%;C:\local\OpenBLAS-v0.2.14-Win32\include
set INCLUDE=%INCLUDE%;C:\local\zlib-1.2.8\include
set INCLUDE=%INCLUDE%;C:\local\MUMPS\include

::set LIB=c:\local\lib


set PATH=c:\MinGW\bin
set PATH=%PATH%;C:\MinGW\msys\1.0\bin
set PATH=%PATH%;C:\local\tbb43_20150424oss\build\windows_ia32_gcc_mingw4.8.1_debug
set PATH=%PATH%;C:\local\tbb43_20150424oss\build\windows_ia32_gcc_mingw4.8.1_release
set PATH=%PATH%;C:\local\gmsh-2.9.3-Windows
set PATH=%PATH%;C:\local\swigwin-2.0.12
set PATH=%PATH%;C:\local\OpenBLAS-v0.2.14-Win32\bin
set PATH=%PATH%;C:\local\qt-4.8.7\bin
set PATH=%PATH%;C:\Python27
set PATH=%PATH%;C:\local\vtk-5.10.1\bin
set PATH=%PATH%;C:\local\zlib-1.2.8\bin
set PATH=%PATH%;c:\local\bin
set PATH=%PATH%;c:\local\MUMPS\bin
set PATH=%PATH%;C:\Program Files\CMake\bin
set PATH=%PATH%;C:\Program Files\TortoiseSVN\bin

%comspec%
