cmake_minimum_required(VERSION 2.6)

project(FEM)

find_package(Armadillo REQUIRED)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -Wall -fcolor-diagnostics")

set( MY_INCLUDES ${PROJECT_SOURCE_DIR}/include)

set( MY_SOURCES src/main.cpp src/Grid.cpp)

include_directories( ${MY_INCLUDES} ${ARMADILLO_INCLUDE_DIRS} )
add_executable(FEM ${MY_SOURCES})
target_link_libraries(FEM ${ARMADILLO_LIBRARIES})
