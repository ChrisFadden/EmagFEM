#ifndef PETSCWRAPPERS_H
#define PETSCWRAPPERS_H

#include <petsc.h>
#include <vector>

Vec *createVec(const std::vector<double> &);

#endif
