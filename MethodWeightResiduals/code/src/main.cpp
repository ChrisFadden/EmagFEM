#include "petscsys.h"
#include "wrapMat.h"
#include "wrapVec.h"
#include <iostream>
#include <petsc.h>
#include <petscviewerhdf5.h>
#include <string>

int main(int argc, char **args) {

  PetscInitialize(&argc, &args, NULL, NULL);
  { // Wrap everything in a scope so destructors are called before PetscFinalize

    // Prevent .info file from being created
    PetscErrorCode ierr =
        PetscOptionsSetValue(NULL, "-viewer_binary_skip_info", "");
    CHKERRQ(ierr);

    std::vector<double> ab = {7.0, 1.0, 1.0, 3.0}; // entries of vector b

    std::vector<std::vector<double>> aA = {
        {1.0, 2.0, 3.0, 0.0}, // entries of matrix A
        {2.0, 1.0, -2.0, -3.0},
        {-1.0, 1.0, 1.0, 0.0},
        {0.0, 1.0, 1.0, -1.0}};

    wrapVec b(ab); // initialize b to vector
    wrapMat A(aA); // initialize A to matrix
    wrapVec x(b);  // initialize answer same as b
    A.solve(b, x);
    b.writeToBIN("vec.dat");
    A.writeToBIN("mat.dat");
  }
  PetscFinalize();
  return 0;
}
