#include "wrapMat.h"
#include <iostream>
#include <petscsys.h>

/***************
 * CONSTRUCTORS
 **************/
// Construct from a std::vector
wrapMat::wrapMat(std::vector<std::vector<real>> &G) {
  PetscErrorCode ierr;

  this->M = G.size();
  this->N = G[0].size();

  ierr = MatCreate(PETSC_COMM_WORLD, &(this->A));
  ierr = MatSetSizes(this->A, PETSC_DECIDE, PETSC_DECIDE, this->M, this->N);
  ierr = MatSetFromOptions(this->A);
  ierr = MatSetUp(this->A);

  std::vector<int> idx(this->N);
  for (int nn = 0; nn != this->N; ++nn) {
    idx[nn] = nn;
  }

  // Set Matrix entries
  for (int mm = 0; mm != this->M; ++mm) {
    ierr = MatSetValues(this->A, 1, &mm, this->N, idx.data(), G[mm].data(),
                        INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(this->A, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(this->A, MAT_FINAL_ASSEMBLY);
}

// Copy Constructor (modify for Sparse Matrices)
wrapMat::wrapMat(wrapMat &old) {
  this->M = old.getSize_M();
  this->N = old.getSize_N();
  MatDuplicate(old.getMat(), MAT_COPY_VALUES, &(this->A));
}

/*********************
 * MEMBER FUNC
 ********************/

// Access Functions
Mat wrapMat::getMat() { return this->A; }
loop wrapMat::getSize_M() { return this->M; }
loop wrapMat::getSize_N() { return this->N; }

//  Matrix Solve
void wrapMat::solve(wrapVec &b, wrapVec &x) {
  PetscErrorCode ierr;
  KSP ksp;
  // Create Solver
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  ierr = KSPSetOperators(ksp, this->A, this->A);
  ierr = KSPSetFromOptions(ksp);

  //  x = b;
  ierr = KSPSolve(ksp, b.getVec(), x.getVec());
  // KSPDestroy(&ksp);

  return;
}

//  Binary File Write
void wrapMat::writeToBIN(std::string file, PetscFileMode type, bool fullPath) {
  // Create the FileName
  std::string fp = "";
  if (!fullPath)
    fp += DATA_PATH;
  fp += file;

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, fp.c_str(), type, &viewer);
  MatView(this->A, viewer);
  PetscViewerDestroy(&viewer);
}

/********************
 * DESTRUCTOR
 *******************/
void wrapMat::cleanMem() { MatDestroy(&(this->A)); }
wrapMat::~wrapMat() { MatDestroy(&(this->A)); }
