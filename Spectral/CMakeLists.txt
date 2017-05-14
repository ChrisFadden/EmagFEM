cmake_minimum_required(VERSION 2.6)

project(SpectralFEM)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -Wall -fcolor-diagnostics")

set( MY_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set( MY_SOURCES src/main.cpp)

include_directories(${MY_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS})
add_executable(SpectralFEM ${MY_SOURCES})

target_link_libraries(SpectralFEM ${OpenCV_LIBS})
