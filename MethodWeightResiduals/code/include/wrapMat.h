#ifndef WRAPMAT_H
#define WRAPMAT_H

#include "Constants.h"
#include "wrapVec.h"
#include <petsc.h>
#include <vector>

class wrapMat {
private:
  Mat A;
  loop M;
  loop N;

public:
  // Constructors/Destructor
  wrapMat(std::vector<std::vector<real>> &);
  wrapMat(wrapMat &);
  void cleanMem();
  ~wrapMat();

  // Access Functions
  Mat getMat();
  loop getSize_M();
  loop getSize_N();

  // Solve
  void solve(wrapVec &b, wrapVec &x);

  void writeToBIN(std::string, PetscFileMode type = FILE_MODE_WRITE,
                  bool fullPath = false);
};

#endif
