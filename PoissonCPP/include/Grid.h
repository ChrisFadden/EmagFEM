#ifndef GRID_H
#define GRID_H

#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <armadillo>

// Structure containing nodes of a square
struct Triangle {
  int a;
  int b;
  int c;
};

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
  std::vector<std::pair<int, int>> nodePairs;
  std::vector<Triangle> elements;
  std::vector<int> boundaryNodes; 
 private:
  void meshNodes();
  void meshElements();
  void setBoundaryNodes();
};

#endif
