#ifndef MATRIXIO_H
#define MATRIXIO_H
#include "Constants.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>

void mat2bin(std::string fname, const Mat &A) {
  std::ofstream fp;
  fp.open(fname, std::ios::out | std::ios::binary);

  auto rows = A.rows();
  auto cols = A.cols();
  auto szByte = sizeof(real);
  fp.write((char *)(&rows), sizeof(rows));
  fp.write((char *)(&cols), sizeof(cols));
  fp.write((char *)(&szByte), sizeof(szByte));
  fp.write((char *)A.data(), rows * cols * sizeof(A(0, 0)));
  fp.close();

  return;
}

void vec2bin(std::string fname, const Vec &x) {
  std::ofstream fp;
  fp.open(fname, std::ios::out | std::ios::binary);

  auto sz = x.size();
  auto szByte = sizeof(real);

  fp.write((char *)(&sz), sizeof(sz));
  fp.write((char *)(&szByte), sizeof(szByte));
  fp.write((char *)x.data(), sz * sizeof(x(0)));
  fp.close();

  return;
}

#endif
