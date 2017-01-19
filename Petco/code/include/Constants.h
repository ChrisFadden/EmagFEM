#ifndef CONSTANTS_H
#define CONSTANTS_H

#define DATA_PATH "../../data/"

#include <complex>

// Precision Numbers (PETSc can only handle double as its compiled (I think...)
typedef double real;
typedef std::complex<real> cplx;

// loop
typedef std::size_t loop;
#endif