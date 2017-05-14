#include "fem1d.h"
#include "linearAlgebraOCV.h"
#include <iostream>

int main(int argc, char **argv) {

  //  testLinearAlgebra();
  // test_fem1d_Diff();

  // test_fem1d_Helm();
  test_fem1d_QuadElm_Diff();

  std::cout << "Hello World" << std::endl;
  return 0;
}
