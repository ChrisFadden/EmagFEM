// Aramdillo headers
#include <armadillo>

// System Headers
#include <iostream>

// Chris typedefs
typedef float real;
typedef unsigned int loop;
typedef unsigned int node;

int main(int argc, char **argv) {
  // Number of points in each direction (M = x, N = y)
  int M = 3;
  int N = 3;
  int numElements = 2 * (M - 1) * (N - 1);

  // Grid Boundaries
  int x0 = 0;
  int x1 = 1;
  real dx = (real)(x1 - x0) / (M - 1);
  int y0 = 0;
  int y1 = 1;
  real dy = (real)(y1 - y0) / (N - 1);

  /********************
   * Triangulate Nodes
   *******************/
  arma::Mat<node> T(numElements, 3);
  // Populate Triangle node numbers
  loop kk = 0;
  for (loop ii = 0; ii != numElements; ++ii) {
    kk = ii / (2 * (M - 1));
    if (ii % 2 == 0) {
      T(ii, 0) = ii / 2 + M + kk;
      T(ii, 1) = ii / 2 + kk;
      T(ii, 2) = ii / 2 + kk + 1;
    } else {
      T(ii, 0) = ii / 2 + M + kk;
      T(ii, 1) = ii / 2 + kk + 1;
      T(ii, 2) = ii / 2 + M + kk + 1;
    }
  } // end element loop

  arma::Mat<real> Tpoints(M * N, 2);
  kk = 0;
  for (loop jj = 0; jj < N; ++jj) {
    for (loop ii = 0; ii < M; ++ii) {
      Tpoints(kk, 0) = ii * dx;
      Tpoints(kk, 1) = jj * dx;
      ++kk;
    }
  }

  /*****************
   * Boundary
   ****************/
  // Boundary nodes
  int boundarySize = 2 * (M - 1) + 2 * (N - 1);
  arma::Col<node> b(boundarySize);
  kk = 0;
  // Bottom
  for (loop ii = 0; ii < M; ++ii) {
    b(kk) = ii;
    ++kk;
  }

  // Left
  for (loop ii = M; ii < M * N; ii += M) {
    b(kk) = ii;
    ++kk;
  }

  // Right
  for (loop ii = 2 * M - 1; ii < M * N; ii += M) {
    b(kk) = ii;
    ++kk;
  }

  // Top
  for (loop ii = M * N - M + 1; ii < M * N - 1; ++ii) {
    b(kk) = ii;
    ++kk;
  }

  /********************
   Form Local Matrices
  ********************/
  arma::Mat<real> Clocal(numElements, 3, 3);
  arma::Col<node> P(3);
  arma::Col<node> Q(3);
  for (loop ii = 0; ii < numElements; ++ii) {
    Tlocal(T(ii, 0), 0); // x of 1st node
    Tlocal(T(ii, 1), 0); // x of 2nd node
    Tlocal(T(ii, 2), 0); // x of 3rd node
    Tlocal(T(ii, 0), 1); // y of 1st node
    Tlocal(T(ii, 1), 1); // y of 2nd node
    Tlocal(T(ii, 2), 1); // y of 3rd node
  }

  // P=zeros(no_of_elems,3);
  // Q=zeros(no_of_elems,3);
  // C=zeros(no_of_elems,3,3);
  // for i1=1:1:no_of_elems
  // P(i1,1)=grid_pts(look_up_table(i1,3),3)-grid_pts(look_up_table(i1,4),3);
  // P(i1,2)=grid_pts(look_up_table(i1,4),3)-grid_pts(look_up_table(i1,2),3);
  // P(i1,3)=grid_pts(look_up_table(i1,2),3)-grid_pts(look_up_table(i1,3),3);
  // Q(i1,1)=grid_pts(look_up_table(i1,4),2)-grid_pts(look_up_table(i1,3),2);
  // Q(i1,2)=grid_pts(look_up_table(i1,2),2)-grid_pts(look_up_table(i1,4),2);
  // Q(i1,3)=grid_pts(look_up_table(i1,3),2)-grid_pts(look_up_table(i1,2),2);
  // for j=1:1:3
  // for k=1:1:3
  // C(i1,j,k)=(1/(4*area_of_elem))*(P(i1,j)*P(i1,k)+Q(i1,j)*Q(i1,k));
  // end
  // end
  // end

  //% Initializing matrix C_global from which propagation voltage distribution
  //% is calculated using iterative analysis
  // C_global=zeros(total_no_pts,total_no_pts);

  //% Calculating matrix C1_global by adding up local C matrices
  // for i=1:1:no_of_elems
  // for j=1:1:3
  // for k=1:1:3
  // C_global(look_up_table(i,j+1),look_up_table(i,k+1))=C_global(look_up_table(i,j+1),look_up_table(i,k+1))+C(i,j,k);
  // end
  // end
  // end

  // std::cout << Tpoints << std::endl;
  return EXIT_SUCCESS;
}
