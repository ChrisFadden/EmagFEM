#ifndef GRID_H
#define GRID_H

#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include <armadillo>

class Grid {
 public:
  // Constructors
  Grid();
  Grid(std::string);  // Filename

  // Functions
  void Mesh();

  // Members
  int sizeX;
  int sizeY;
  int numNodes;
  std::vector<std::pair<double, double>> nodePairs;
  arma::mat triangles;
};

#endif
