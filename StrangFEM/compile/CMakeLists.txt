cmake_minimum_required(VERSION 2.6)
project(FEM)

find_package (Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(CGAL QUIET COMPONENTS Core )
include(${CGAL_USE_FILE})

#Add C++ Flags
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -Wall -fcolor-diagnostics")

#Set Sources
set(CODE_DIR ${CMAKE_SOURCE_DIR}/../)
include_directories(${CODE_DIR}/include) 

#Create and Link Executable
add_executable(FEM ${CODE_DIR}/src/main.cpp)
target_link_libraries(FEM ${ARMADILLO_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
