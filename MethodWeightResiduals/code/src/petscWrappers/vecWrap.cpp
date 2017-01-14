#include "petscWrappers.h"
#include <iostream>
#include <stdlib.h>

typedef unsigned loop;

// Wrapper for creating a PETSc vector from a C++ vector
Vec *createVec(const std::vector<double> &v) {
  PetscErrorCode ierr;

  Vec *b = (Vec *)malloc(sizeof(Vec));
  std::vector<int> idx(v.size());
  for (loop i = 0; i != idx.size(); ++i) {
    idx[i] = i;
  }

  ierr = VecCreate(PETSC_COMM_WORLD, b);
  ierr = VecSetSizes(*b, PETSC_DECIDE, v.size());
  ierr = VecSetFromOptions(*b);
  ierr = VecSetValues(*b, v.size(), idx.data(), v.data(), INSERT_VALUES);

  // Distribute with MPI
  ierr = VecAssemblyBegin(*b);
  ierr = VecAssemblyEnd(*b);
  return b;
}
