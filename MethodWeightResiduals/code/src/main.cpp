#include "petscsys.h"
#include "wrapMat.h"
#include "wrapVec.h"
#include <iostream>
#include <petsc.h>
#include <petscviewerhdf5.h>

int main(int argc, char **args) {
  PetscErrorCode ierr;
  PetscInitialize(&argc, &args, NULL, NULL);
  { // Wrap everything in a scope so destructors are called before PetscFinalize
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

    // Print Out Vector
    ierr = VecView(x.getVec(), PETSC_VIEWER_STDOUT_WORLD);

    // Print Vec to file
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, "../../data/vector.h5",
                        FILE_MODE_WRITE, &viewer);
    VecView(x.getVec(), viewer);
    PetscViewerDestroy(&viewer);
  }
  PetscFinalize();
  return 0;
}
