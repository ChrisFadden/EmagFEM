#ifndef WRAPVEC_H
#define WRAPVEC_H

#include "Constants.h"
#include <petsc.h>
#include <petscviewerhdf5.h>
#include <string>
#include <vector>
class wrapVec {
private:
  Vec v;
  loop sizeVec;

public:
  // Constructors/Destructor
  wrapVec(std::vector<real> &);
  wrapVec(wrapVec &);
  ~wrapVec();
  void cleanMem();

  // Access Functions
  Vec getVec();
  loop getSize();

  // Print Functions
  void writeToBIN(std::string, PetscFileMode type = FILE_MODE_WRITE,
                  bool fullPath = false);
};

#endif
