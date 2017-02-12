#include <iostream>
#include <vector>

#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/element_creation.hpp"
#include "viennagrid/mesh/mesh.hpp"
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

  // Construct mesh object
  MeshType mesh;

  // Define points
  VertexHandleType vh0 =
      viennagrid::make_vertex(mesh, PointType(0.8, 1.8)); // id = 0
  VertexHandleType vh1 =
      viennagrid::make_vertex(mesh, PointType(1.4, 1.4)); // id = 1
  VertexHandleType vh2 = viennagrid::make_vertex(mesh, PointType(2.1, 2.1));
  VertexHandleType vh3 = viennagrid::make_vertex(mesh, PointType(1.2, 2.7));

  // Create triangles
  viennagrid::make_triangle(mesh, vh0, vh1, vh3);
  viennagrid::make_triangle(mesh, vh1, vh2, vh3);

  // Populate Stiffness matrix
  for (auto elem = viennagrid::elements<viennagrid::triangle_tag>(mesh).begin();
       elem != viennagrid::elements<viennagrid::triangle_tag>(mesh).end();
       elem++) {

    auto vert = viennagrid::elements<viennagrid::vertex_tag>(*elem);
    auto v0 = viennagrid::point(vert[0]);
    auto v1 = viennagrid::point(vert[1]);
    auto v2 = viennagrid::point(vert[2]);

    float P[3];
    float Q[3];

    P[0] = (v1[1] - v2[1]);
    P[1] = (v2[1] - v0[1]);
    P[2] = (v0[1] - v1[1]);

    Q[0] = (v2[0] - v1[0]);
    Q[1] = (v0[0] - v2[0]);
    Q[2] = (v1[0] - v0[0]);

    float A = 0.5 * (P[1] * Q[2] - P[2] * Q[1]);
    float Cfact = 1 / (4 * A);
    for (int ii = 0; ii < 3; ++ii) {
      for (int jj = 0; jj < 3; ++jj) {
        std::cout << Cfact * (P[ii] * P[jj] + Q[ii] * Q[jj]) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
