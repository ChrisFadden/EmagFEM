#ifndef PETCO_H
#define PETCO_H

#include "Constants.h"
#include <petsc.h>
#include <string>
#include <vector>

static PetscErrorCode ierr;

namespace petco {

// Construct Vec from 1D std::vector
void VecFromVector(std::vector<real> &, Vec *);

// Write Vector to binary file
void writeVecToBIN(Vec, std::string, PetscFileMode type = FILE_MODE_WRITE,
                   bool fullPath = false);

// Construct Mat from 2D std::vector<std::vector<real>>
void MatFromVector(std::vector<std::vector<real>> &, Mat *);

// Write Matrix to binary file
void writeMatToBIN(Mat, std::string, PetscFileMode type = FILE_MODE_WRITE,
                   bool fullPath = false);
}

#endif
