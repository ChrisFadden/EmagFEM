#include <iostream>
#include <armadillo>
#include <cmath>

#include "Grid.h"

void applyBC(arma::mat,arma::vec,Grid&);

int main() {
  Grid g;
  g.Mesh();

  // Execute Finite Elements
  arma::mat K(g.nodePairs.size(), g.nodePairs.size());
  arma::vec F(g.nodePairs.size());

  for (int i = 0; i < g.triangles.n_rows; i++) {
    arma::mat elem;
    elem << 1 << g.nodePairs[g.triangles(i, 0)].first << g.nodePairs[g.triangles(i, 0)].second << arma::endr 
       << 1 << g.nodePairs[g.triangles(i, 1)].first << g.nodePairs[g.triangles(i, 1)].second << arma::endr 
       << 1 << g.nodePairs[g.triangles(i, 2)].first << g.nodePairs[g.triangles(i, 2)].second << arma::endr;
  
    double area = std::abs(arma::det(elem)) / 2.0; 
    
    arma::mat C = arma::inv(elem);
    arma::mat grad = C.rows(1,2);
    
    arma::mat Ke = area*trans(grad)*grad;
    double Fe = area / 3.0;

    for(int j = 0; j < 3; j++)
    {
      for(int k = 0; k < 3; k++)
      {
        K(g.triangles(i,j),g.triangles(i,k)) += Ke(j,k); 
      }
      F(g.triangles(i,j)) += Fe;
    }
  }
  
  applyBC(K,F,g);

  arma::vec U = solve(K,F);
  
  std::cout << U << std::endl;

  std::cout << "Hello World" << std::endl;

  return 0;
}

//set Dirichlet BC
void applyBC(arma::mat K, arma::vec F, Grid& g){
  std::vector<int> b(2*(g.sizeX-1) + 2*(g.sizeY-1));
  
  int i = 0;
  for(int j = 0; j < g.sizeX; j++)
  {
    b[i] = j;
    i++;
  }
  for(int j = g.sizeX; j < g.sizeX*g.sizeY; j+=g.sizeX){
    b[i] = j;
    i++; 
  }

  for(int j = 2*g.sizeX-1; j < g.sizeX*g.sizeY; j+=g.sizeX){
    b[i] = j;
    i++;
    std::cout << j << std::endl;
  }
  
  for(int j = g.sizeX*g.sizeY - g.sizeX + 1; j < g.sizeX*g.sizeY - 1; j++){
    b[i] = j;
    i++;
  }
  
  for(auto i : b){ 
    K.row(i) = 0 * K.row(i);
    K.col(i) = 0 * K.col(i); 
    F.at(i) = 0; 
  }
}



