cmake_minimum_required(VERSION 2.6)

project(FEM)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -Wall -fcolor-diagnostics")

set( MY_INCLUDES ${PROJECT_SOURCE_DIR}/include)

set( MY_SOURCES src/main.cpp)

include_directories( ${MY_INCLUDES} /usr/include/eigen3 )
add_executable(FEM ${MY_SOURCES})
