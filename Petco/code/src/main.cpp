#include "petco.h"
#include "petscsys.h"
#include <iostream>
#include <petsc.h>
#include <petscviewerhdf5.h>
#include <string>

int main(int argc, char **args) {

  PetscInitialize(&argc, &args, NULL, NULL);
  // Prevent .info file from being created
  ierr = PetscOptionsSetValue(NULL, "-viewer_binary_skip_info", "");
  CHKERRQ(ierr);

  std::vector<double> ab = {7.0, 1.0, 1.0, 3.0}; // entries of vector b

  std::vector<std::vector<double>> aA = {
      {1.0, 2.0, 3.0, 0.0}, // entries of matrix A
      {2.0, 1.0, -2.0, -3.0},
      {-1.0, 1.0, 1.0, 0.0},
      {0.0, 1.0, 1.0, -1.0}};

  Vec x, b;
  Mat A;
  KSP ksp;

  petco::VecFromVector(ab, &b);
  petco::MatFromVector(aA, &A);

  // Solve Matrix
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  ierr = KSPSetOperators(ksp, A, A);
  ierr = KSPSetFromOptions(ksp);
  ierr = VecDuplicate(b, &x);
  ierr = KSPSolve(ksp, b, x);
  petco::writeVecToBIN(x, "vec.dat");
  petco::writeMatToBIN(A, "mat.dat");

  KSPDestroy(&ksp);
  VecDestroy(&b);
  VecDestroy(&x);
  MatDestroy(&A);

  PetscFinalize();
  return 0;
}
