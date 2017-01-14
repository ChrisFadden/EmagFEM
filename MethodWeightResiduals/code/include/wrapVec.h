#ifndef WRAPVEC_H
#define WRAPVEC_H

#include "Constants.h"
#include <petsc.h>
#include <vector>

class wrapVec {
private:
  Vec v;
  loop sizeVec;

public:
  // Constructors/Destructor
  wrapVec(const std::vector<real> &);
  void cleanMem();

  // Access Functions
  Vec getVec();
};

#endif
