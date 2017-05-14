#include <armadillo>
#include <iostream>

int main() {
  arma::vec b;
  b << 2.0 << 5.0 << 2.0;

  // arma::endr represents the end of a row in a matrix
  arma::mat A;
  A << 1.0 << 2.0 << arma::endr << 2.0 << 3.0 << arma::endr << 1.0 << 3.0
    << arma::endr;

  std::cout << "Least squares solution:\n";
  std::cout << arma::solve(A, b) << '\n';

  return 0;
}
