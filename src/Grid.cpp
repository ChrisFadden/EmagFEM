#include "Grid.h"
#include <armadillo>

Grid::Grid() {
  sizeX = 5;
  sizeY = 5;
}

Grid::Grid(std::string filename) {
  sizeX = 3;
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
  elements.resize(2*(sizeX-1)*(sizeY-1)); 
  int j = 0;
  int k = 0;
  int l = 0;
  for(int i = 0; i < elements.size(); i++){
    elements[i].a = (sizeX)*j + k;
    elements[i].b = (sizeX)*(j+l) + 1 + k;
    elements[i].c = (sizeX)*(j+1) + 1 + k - l;
    j++; 
    if(j / float(sizeY-1) >= 1) 
    {
      j = 0;
      k++;
    }
    
    if(k == sizeX-1){
      k = 0;
      j = 0;
      l++;
    }

    std::cout << elements[i].a+1 << " " << elements[i].b+1 << " " << elements[i].c + 1 << std::endl;
  }
  
}  // end mesh Elements
