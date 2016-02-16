#include "Grid.h"
#include <iostream>
#include <armadillo>

Grid::Grid() {
  sizeX = 4;
  sizeY = 4;
}

Grid::Grid(std::string filename) {
  sizeX = 4;
  sizeY = 4;
  std::string notUsed = filename;
}

void Grid::Mesh() {
  // Create Node Pairs
  double x, y;
  nodePairs.resize(sizeX * sizeY);

  int node = 0;

  for (int j = 0; j < sizeY; j++) {
    y = (double)j / (sizeY);
    for (int i = 0; i < sizeX; i++) {
      x = (double)i / (sizeX);
      nodePairs[node].first = x;
      nodePairs[node].second = y;
      node++;
    }  // end sizeY
  }  // end sizeX

  // Create Triangles (Check if I subtracted 1 correctly...)
  arma::mat tmpMatrix1, tmpMatrix3;
  arma::vec tmpVec1, tmpVec2;

  // X vertices
  tmpMatrix1 << 1 << 2 << sizeX + 2 << arma::endr << 1 << sizeX + 2 << sizeX+1
             << arma::endr;
  tmpVec1.ones(sizeX - 1);

  arma::mat tmpMatrix2(size(tmpMatrix1), arma::fill::ones);
  tmpVec2 = arma::linspace<arma::vec>(0, sizeX - 2, sizeX - 1);

  tmpMatrix3 =
      arma::kron(tmpMatrix1, tmpVec1) + arma::kron(tmpMatrix2, tmpVec2);

  // Y vertices
  arma::vec tmpVec3, tmpVec4;
  tmpVec3.ones(sizeY - 1);

  arma::mat tmpMatrix4(size(tmpMatrix3), arma::fill::ones);
  tmpVec4 = arma::linspace<arma::vec>(0, sizeY - 2, sizeY - 1);
  tmpVec4 = tmpVec4 * sizeX;

  triangles = arma::kron(tmpMatrix3, tmpVec3) + arma::kron(tmpMatrix4, tmpVec4);
  triangles = triangles - 1;
}  // end Mesh func
