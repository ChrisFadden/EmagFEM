cmake_minimum_required(VERSION 2.6)
project(MWR)

#Find PETSC
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")
find_package (PETSc REQUIRED)
include_directories (${PETSC_INCLUDES})
add_definitions (${PETSC_DEFINITIONS})

#Add C++ Flags
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -Wall -fcolor-diagnostics")

#Set Sources
set(CODE_DIR ${CMAKE_SOURCE_DIR}/../code)
include_directories(${CODE_DIR}/include) 

add_subdirectory(${CODE_DIR}/src/petco ${CMAKE_SOURCE_DIR}/build/petco)

#Create and Link Executable
add_executable(MWR ${CODE_DIR}/src/main.cpp)
target_link_libraries(MWR ${PETSC_LIBRARIES} petcoLib)
