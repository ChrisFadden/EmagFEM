#include <iostream>
#include <armadillo>
#include <cmath>

#include "Grid.h"

int main() {
  Grid g;
  g.Mesh();
  
  //Form Finite Element matrices
  int N = g.sizeX*g.sizeY; 
  arma::mat K(N,N,arma::fill::zeros);
  arma::vec F(N,arma::fill::zeros);
  
  for(int i = 0; i < g.elements.size(); i++){
    Triangle nodes = g.elements[i];
    arma::mat Pe(3,3,arma::fill::ones);
    
    for(int i = 0; i < 3; i++){
      int node; 
      if(i == 0){
        node = nodes.a;
      } else if (i == 1){
        node = nodes.b;
      } else{
        node = nodes.c;
      }
      
      Pe(i,1) = g.nodePairs[node].first;
      Pe(i,2) = g.nodePairs[node].second;
    }
    
    //Find area of the triangle elements
    double Area = std::abs(arma::det(Pe)) / 2.0;
    
    arma::mat C = arma::inv(Pe); 
    arma::mat grad = C.rows(1,2);
    
    arma::mat Ke = Area * grad.t() * grad;
    double Fe = Area / 3.0;

    //Update Global Matrices
    K(nodes.a,nodes.a) += Ke(0,0);
    K(nodes.a,nodes.b) += Ke(0,1);
    K(nodes.a,nodes.c) += Ke(0,2);

    K(nodes.b,nodes.a) += Ke(1,0);
    K(nodes.b,nodes.b) += Ke(1,1);
    K(nodes.b,nodes.c) += Ke(1,2);
    
    K(nodes.c,nodes.a) += Ke(2,0);
    K(nodes.c,nodes.b) += Ke(2,1);
    K(nodes.c,nodes.c) += Ke(2,2);

    F(nodes.a) += Fe;
    F(nodes.b) += Fe;
    F(nodes.c) += Fe;
  }
  
  //Apply Boundary Conditions
  for(int i = 0; i < g.boundaryNodes.size(); i++)
    K.row(g.boundaryNodes[i]) = arma::rowvec(N,arma::fill::zeros);
  for(int i = 0; i < g.boundaryNodes.size(); i++)
    K.col(g.boundaryNodes[i]) = arma::vec(N,arma::fill::zeros);
  for(int i = 0; i < g.boundaryNodes.size(); i++)
    K(g.boundaryNodes[i],g.boundaryNodes[i]) = 1;
  for(int i = 0; i < g.boundaryNodes.size();i++)
    F(g.boundaryNodes[i]) = 0;

  //Find Solution
  arma::vec U = arma::solve(K,F);
  
  std::cout << U << std::endl;

  return 0;
}

