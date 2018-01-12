cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(NRG_seminar)

find_package(PCL 1.7.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (oseg ObjectSegmentation.hpp ObjectSegmentation.cpp CSCSegmentation.cpp)
target_link_libraries (oseg ${PCL_LIBRARIES})
