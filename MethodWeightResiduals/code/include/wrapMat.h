#ifndef WRAPMAT_H
#define WRAPMAT_H

#include "Constants.h"
#include <petsc.h>
#include <vector>

class wrapMat {
private:
  Mat A;
  loop M;
  loop N;

public:
  // Constructors/Destructor
  wrapMat(const std::vector<std::vector<real>> &);
  wrapMat(wrapMat &);
  void cleanMem();

  // Access Functions
  Mat getMat();
  loop getSize_M();
  loop getSize_N();
};

#endif
