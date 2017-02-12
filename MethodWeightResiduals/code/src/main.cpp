// STL
#include <iostream>
#include <vector>

// GRID
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/element_creation.hpp"
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/point.hpp"

// Linear Algebra
#include "viennacl/linalg/direct_solve.hpp"
#include "viennacl/linalg/lu.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/vector.hpp"
#include "viennagrid/point.hpp"

typedef viennagrid::triangular_2d_mesh MeshType;
typedef viennagrid::result_of::segmentation<MeshType>::type SegmentationType;
typedef viennagrid::result_of::segment_handle<SegmentationType>::type
    SegmentHandleType;

typedef viennagrid::result_of::point<MeshType>::type PointType;
typedef viennagrid::result_of::vertex_handle<MeshType>::type VertexHandleType;

typedef viennagrid::result_of::cell<MeshType>::type CellType;
typedef viennagrid::result_of::cell_range<SegmentHandleType>::type CellRange;
typedef viennagrid::result_of::iterator<CellRange>::type CellIterator;

int main() {

  /***************
   * CONSTANTS
   **************/
  MeshType mesh;
  const int numX = 2;
  const int numY = 2;
  int numElements = 2 * numX * numY;
  int numPoints = (numX + 1) * (numY + 1);

  /******************
   * MESHING
   *****************/
  // populate the triangles
  VertexHandleType vh[numPoints];
  int kk = 0;
  for (int jj = 0; jj < numY + 1; ++jj) {
    for (int ii = 0; ii < numX + 1; ++ii) {
      vh[kk] = viennagrid::make_vertex(mesh, PointType(ii, jj));
      ++kk;
    }
  }

  for (int ii = 0; ii != numElements; ++ii) {
    kk = ii / (2 * numX);
    if (ii % 2 == 0) {
      viennagrid::make_triangle(mesh, vh[ii / 2 + numX + kk + 1],
                                vh[ii / 2 + kk], vh[ii / 2 + kk + 1]);
    } else {
      viennagrid::make_triangle(mesh, vh[ii / 2 + numX + kk + 1],
                                vh[ii / 2 + kk + 1],
                                vh[ii / 2 + numX + kk + 2]);
    }
  } // end element loop

  /******************
   * Stiffness Matrix
   *****************/
  float CLocal[numElements][3][3];
  float P[3];
  float Q[3];
  float A = 0.5; // area of element (1/2 base * height)
  for (int ii = 0; ii != numElements; ++ii) {
    auto elem = viennagrid::elements<viennagrid::triangle_tag>(mesh)[ii];
    auto v1 = viennagrid::point(
        viennagrid::elements<viennagrid::vertex_tag>(elem)[0]);
    auto v2 = viennagrid::point(
        viennagrid::elements<viennagrid::vertex_tag>(elem)[1]);
    auto v3 = viennagrid::point(
        viennagrid::elements<viennagrid::vertex_tag>(elem)[2]);

    P[0] = v2[1] - v3[1];
    P[1] = v3[1] - v1[1];
    P[2] = v1[1] - v2[1];
    Q[0] = v3[0] - v2[0];
    Q[1] = v1[0] - v3[0];
    Q[2] = v2[0] - v1[0];

    for (int jj = 0; jj != 3; ++jj) {
      for (int kk = 0; kk != 3; ++kk) {
        CLocal[ii][jj][kk] = (1 / (4 * A)) * (P[jj] * P[kk] + Q[jj] * Q[kk]);
      }
    }
  }

  viennacl::matrix<float> CGlobal(numPoints, numPoints);
  auto elem = viennagrid::elements<viennagrid::triangle_tag>(mesh)[0];
  int vertX;
  int vertY;
  for (int ii = 0; ii != numElements; ++ii) {
    elem = viennagrid::elements<viennagrid::triangle_tag>(mesh)[ii];
    for (int jj = 0; jj != 3; ++jj) {
      vertX = viennagrid::elements<viennagrid::vertex_tag>(elem)[jj].id().get();
      for (int kk = 0; kk != 3; ++kk) {
        vertY =
            viennagrid::elements<viennagrid::vertex_tag>(elem)[kk].id().get();
        CGlobal(vertX, vertY) += CLocal[ii][jj][kk];
      }
    }
  }

  std::cout << CGlobal << std::endl;

  return EXIT_SUCCESS;
}
