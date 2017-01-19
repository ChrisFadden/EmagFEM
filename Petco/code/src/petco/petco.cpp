#include "petco.h"
#include <iostream>

/***************
 * CONSTRUCTORS
 **************/

void petco::VecFromVector(std::vector<real> &b, Vec *v) {
  //  PetscErrorCode ierr;

  std::vector<int> idx(b.size());
  for (loop i = 0; i != idx.size(); ++i) {
    idx[i] = i;
  }

  ierr = VecCreate(PETSC_COMM_WORLD, v);
  ierr = VecSetSizes(*v, PETSC_DECIDE, b.size());
  ierr = VecSetFromOptions(*v);
  ierr = VecSetValues(*v, b.size(), idx.data(), b.data(), INSERT_VALUES);

  // Distribute with MPI
  ierr = VecAssemblyBegin(*v);
  ierr = VecAssemblyEnd(*v);
}

void petco::MatFromVector(std::vector<std::vector<real>> &G, Mat *A) {
  // PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, A);
  ierr = MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, G.size(), G[0].size());
  ierr = MatSetFromOptions(*A);
  ierr = MatSetUp(*A);

  std::vector<int> idx(G[0].size());
  for (int nn = 0; nn != G[0].size(); ++nn) {
    idx[nn] = nn;
  }

  // Set Matrix entries
  for (int mm = 0; mm != G.size(); ++mm) {
    ierr = MatSetValues(*A, 1, &mm, G[0].size(), idx.data(), G[mm].data(),
                        INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
}

void petco::writeVecToBIN(Vec v, std::string file, PetscFileMode type,
                          bool fullPath) {
  // Create the FileName
  std::string fp = "";
  if (!fullPath)
    fp += DATA_PATH;
  fp += file;

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, fp.c_str(), type, &viewer);
  VecView(v, viewer);
  PetscViewerDestroy(&viewer);
}

//  Binary File Write
void petco::writeMatToBIN(Mat A, std::string file, PetscFileMode type,
                          bool fullPath) {
  // Create the FileName
  std::string fp = "";
  if (!fullPath)
    fp += DATA_PATH;
  fp += file;

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, fp.c_str(), type, &viewer);
  MatView(A, viewer);
  PetscViewerDestroy(&viewer);
}
