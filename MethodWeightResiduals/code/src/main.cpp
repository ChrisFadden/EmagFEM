#include "wrapVec.h"
#include <iostream>
#include <petsc.h>

int main(int argc, char **args) {
  PetscErrorCode ierr;
  PetscInitialize(&argc, &args, NULL, NULL);

  Vec x;
  Mat A;
  KSP ksp;
  int i, j[4] = {0, 1, 2, 3}; // j = column indices
  std::vector<double> ab = {7.0, 1.0, 1.0, 3.0};

  double aA[4][4] = {{1.0, 2.0, 3.0, 0.0}, // entries of matrix A
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

  // Create Matrix
  ierr = MatCreate(PETSC_COMM_WORLD, &A);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 4, 4);
  ierr = MatSetFromOptions(A);
  ierr = MatSetUp(A);

  { // Set Matrix entries
    for (i = 0; i < 4; i++) {
      ierr = MatSetValues(A, 1, &i, 4, j, aA[i], INSERT_VALUES);
    }
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  }

  /****************
   * MATRIX_SOLVE
   ****************/

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  ierr = KSPSetOperators(ksp, A, A);
  ierr = KSPSetFromOptions(ksp);
  ierr = VecDuplicate(b.getVec(), &x);
  ierr = KSPSolve(ksp, b.getVec(), x);
  ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  /******************
   * CLEAN UP MEMORY
   ****************/

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&x);
  b.cleanMem();

  PetscFinalize();
  return 0;
}
