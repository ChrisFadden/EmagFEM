#include "Grid.h"
#include <iostream>
#include <armadillo>

Grid::Grid() {
  sizeX = 4;
  sizeY = 3;
}

Grid::Grid(std::string filename) {
  sizeX = 4;
  sizeY = 3;
  std::string notUsed = filename;
}

void Grid::Mesh() {
  meshNodes();
  meshElements();
}  // end Mesh func

void Grid::meshNodes() {
  // Create Node Pairs
  nodePairs.resize(sizeX * sizeY);

  int node = 0;

  for (int j = 0; j < sizeY; j++) {
    for (int i = 0; i < sizeX; i++) {
      nodePairs[node].first = i;
      nodePairs[node].second = j;
      node++;
    }  // end sizeY
  }    // end sizeX

}  // end meshNodes function

void Grid::meshElements() {
  // Assumes column major indexing (node # increase along x-axis)
  elements.reserve(sizeX * sizeY);
  Triangle element;

  // Boundary Condition data structures
  std::vector<int> leftBCvec;
  std::vector<int> rightBCvec;
  std::vector<int> downBCvec;
  std::vector<int> upBCvec;

  leftBCvec.reserve(sizeY);
  rightBCvec.reserve(sizeY);
  downBCvec.reserve(sizeX);
  upBCvec.reserve(sizeX);
  boundaryNode.insert({"Left", leftBCvec});
  boundaryNode.insert({"Right", rightBCvec});
  boundaryNode.insert({"Down", downBCvec});
  boundaryNode.insert({"Up", upBCvec});

  int a, b;
  for (int i = 0; i < nodePairs.size(); i++) {
    a = nodePairs[i].first;
    b = nodePairs[i].second;

    if (a > 0 && b + 1 < sizeY) {
      element.a = i;
      element.b = i + sizeX;
      element.c = i - 1 + sizeX;
      elements.push_back(element);
    }

    if (a + 1 < sizeX && b + 1 < sizeY) {
      element.a = i;
      element.b = i + 1;
      element.c = i + sizeX;
      elements.push_back(element);
    }
  }  // end node loop

  // Left BC
  for (int i = 0; i < nodePairs.size(); i += sizeX)
    boundaryNode["Left"].push_back(i);

  // Right BC
  for (int i = sizeX - 1; i < nodePairs.size(); i += sizeX)
    boundaryNode["Right"].push_back(i);

  // Down BC
  for (int i = 0; i < sizeX; i++) boundaryNode["Down"].push_back(i);

  // Up BC
  for (int i = sizeX * (sizeY - 1); i < nodePairs.size(); i++)
    boundaryNode["Up"].push_back(i);

}  // end mesh Elements
