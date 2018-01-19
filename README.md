# Point cloud segmentation using spectral analysis

Point cloud segmentation is implemented in C++. It relies on the following llibraries:
* Eigen
* pcl

## Project structure
* CMakeLists.txt - CMake file
* ObjectSegmentation.hpp - main header file
* ObjectSegmentation.cpp - entry file
* CSCSegmentation.cpp - constrained supervoxel clustering method

## Input
Executable expects files of type .ply
./oseg  pathToFile/pointCloud.ply

## Output
Executable produces files of type .ply
For each segment a file is produced according to label and saved in current directory.