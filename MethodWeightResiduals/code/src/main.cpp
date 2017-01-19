#include <iostream>
#include <vector>

#include "viennacl/linalg/direct_solve.hpp"
#include "viennacl/linalg/lu.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/vector.hpp"

int main() {

  typedef float ScalarType;

  viennacl::vector<ScalarType> b;
  viennacl::vector<ScalarType> x;
  viennacl::matrix<ScalarType> A(4, 4);

  std::vector<double> ab = {7.0, 1.0, 1.0, 3.0}; // entries of vector b

  std::vector<std::vector<double>> aA = {
      {1.0, 2.0, 3.0, 0.0}, // entries of matrix A
      {2.0, 1.0, -2.0, -3.0},
      {-1.0, 1.0, 1.0, 0.0},
      {0.0, 1.0, 1.0, -1.0}};

  viennacl::copy(ab, b);

  for (int i = 0; i != 4; ++i) {
    for (int j = 0; j != 4; ++j) {
      A(i, j) = aA[i][j];
    }
  }

  viennacl::linalg::lu_factorize(A);
  x = b;
  viennacl::linalg::lu_substitute(A, x);

  std::cout << x << std::endl;

  return EXIT_SUCCESS;
}
