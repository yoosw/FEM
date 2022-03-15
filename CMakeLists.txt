cmake_minimum_required(VERSION 3.21)
project(FEM)

set(CMAKE_CXX_STANDARD 14)


add_executable(FEM FEM.cpp ReadInput.cpp loadCrossSection.cpp input-1d.txt CrossSection.txt elements.cpp UpdateKeff.cpp MakeMatrix.cpp)

