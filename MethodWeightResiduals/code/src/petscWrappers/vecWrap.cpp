#include "wrapVec.h"
#include <iostream>
/***************
 * CONSTRUCTORS
 **************/

// Construct from a std::vector
wrapVec::wrapVec(const std::vector<real> &b) {
  PetscErrorCode ierr;

  std::vector<int> idx(b.size());
  for (loop i = 0; i != idx.size(); ++i) {
    idx[i] = i;
  }

  this->sizeVec = b.size();

  ierr = VecCreate(PETSC_COMM_WORLD, &(this->v));
  ierr = VecSetSizes(this->v, PETSC_DECIDE, b.size());
  ierr = VecSetFromOptions(this->v);
  ierr = VecSetValues(this->v, b.size(), idx.data(), b.data(), INSERT_VALUES);

  // Distribute with MPI
  ierr = VecAssemblyBegin(this->v);
  ierr = VecAssemblyEnd(this->v);
}

/*********************
 * MEMBER FUNC
 ********************/
Vec wrapVec::getVec() { return this->v; }

/********************
 * DESTRUCTOR
 *******************/
void wrapVec::cleanMem() { VecDestroy(&v); }
