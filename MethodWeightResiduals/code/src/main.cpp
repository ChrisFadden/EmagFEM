#include "wrapMat.h"
#include "wrapVec.h"
#include <iostream>
#include <petsc.h>

int main(int argc, char **args) {
  PetscErrorCode ierr;
  PetscInitialize(&argc, &args, NULL, NULL);

  KSP ksp;
  std::vector<double> ab = {7.0, 1.0, 1.0, 3.0}; // entries of vector b

  std::vector<std::vector<double>> aA = {
      {1.0, 2.0, 3.0, 0.0}, // entries of matrix A
      {2.0, 1.0, -2.0, -3.0},
      {-1.0, 1.0, 1.0, 0.0},
      {0.0, 1.0, 1.0, -1.0}};

  /******************
   * VECTOR CREATION
   ******************/
  wrapVec b(ab);

  /*****************
   * MATRIX CREATION
   ****************/
  wrapMat A(aA);

  /****************
   * MATRIX_SOLVE
   ****************/

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  ierr = KSPSetOperators(ksp, A.getMat(), A.getMat());
  ierr = KSPSetFromOptions(ksp);

  wrapVec x(b);
  ierr = KSPSolve(ksp, b.getVec(), x.getVec());

  // Print Out Vector
  ierr = VecView(x.getVec(), PETSC_VIEWER_STDOUT_WORLD);

  /******************
   * CLEAN UP MEMORY
   ****************/

  KSPDestroy(&ksp);
  // MatDestroy(&A);
  A.cleanMem();
  x.cleanMem();
  b.cleanMem();

  PetscFinalize();
  return 0;
}
