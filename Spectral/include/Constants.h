#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <Eigen/Dense>
#include <assert.h>
#include <cstddef>
#include <iostream>

#define DBGMODE

#ifdef DBGMODE
#define DEBUG(x)                                                               \
  do {                                                                         \
    std::cout << #x << ":" << std::endl << x << std::endl;                     \
  } while (0);
#else
#define DEBUG(x)
#endif

typedef float real;
typedef std::complex<real> cplx;
typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> Mat;
typedef Eigen::Matrix<real, Eigen::Dynamic, 1> Vec;

typedef std::size_t loop;
typedef std::size_t node;

#endif
