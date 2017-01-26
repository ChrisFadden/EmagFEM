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

  // Print out values
  for (auto iter = viennagrid::elements<viennagrid::triangle_tag>(mesh).begin();
       iter != viennagrid::elements<viennagrid::triangle_tag>(mesh).end();
       iter++) {

    std::cout << *iter << std::endl;
  }

  return EXIT_SUCCESS;
}
