#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/number_utils.h>
#include <math.h>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_ratio_stop_predicate.h>

#include <chrono>
#include <fstream>
#include <iostream>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3>            Surface_mesh;
typedef Surface_mesh::Vertex_index             Vertex_index;
typedef Surface_mesh::Edge_index               Edge_index;
typedef Surface_mesh::Halfedge_index           Halfedge_index;
typedef Surface_mesh::Face_range               Face_range;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor Vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::vertex_iterator   Vertex_iterator;
typedef boost::graph_traits<Surface_mesh>::edge_iterator     Edge_iterator;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor   Edge_descriptor;

//namespace
namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

//height function acting as utility to define the underlying surface
double height(double x, double y) {
  return (sin(x)+sin(y))/2.0;
}

//we can define our surface implicitly via a constraint -- mirroring sphere example
FT heightmap_function(Point_3 p) {
  const FT x = p.x(), y = p.y(), z = p.z();
  return z - height(x,y);
}

float meanEdgeLength(Surface_mesh mesh){
  float mean = 0, min, max, length;
  int count = 0; bool init = true;
  Surface_mesh::Halfedge_range es = mesh.halfedges();
  for (auto  eIter = es.begin(); eIter != es.end(); ++eIter){
      Halfedge_index e = *eIter; 

      Point_3 a = mesh.point(mesh.source(e)); //source yields the source vertex of an edge as a vertex_descriptor
      Point_3 b = mesh.point(mesh.target(e));  

      length = CGAL::sqrt(CGAL::squared_distance(a, b));
      ++count;
      if (init){
          mean = min = max = length;
          init = false;
      }
      else{
          if (length < min) min = length;
          if (length > max) max = length;
      }
      mean += length;
  }
  mean /= count;
  std::cout << min << " " << max << " " << mean << "\n";
  return mean;
}

int main() {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
  // defining the surface
  Surface_3 surface(heightmap_function,             // pointer to function
                    Sphere_3(CGAL::ORIGIN, 20.)); // bounding sphere
  // Note that "20." above is the *squared* radius of the bounding sphere, which
  // I believe defines the domain of the surface
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     0.1,  // radius bound
                                                     0.1); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
  Surface_mesh sm;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
  std::ofstream out("heightmap.off");
  out << sm << std::endl;
  std::cout << "Final number of points pre-simplification: " << tr.number_of_vertices() << "\n";

  //simplify mesh
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  // In this example, the simplification stops when the number of undirected edges
  // drops below 75% of the initial count
  double stop_ratio = 0.60;
  SMS::Edge_count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);
  int r = SMS::edge_collapse(sm, stop);
  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
  std::cout << "\nFinished!\n" << r << " edges removed.\n" << sm.number_of_edges() << " final edges.\n";
  std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;
  CGAL::IO::write_polygon_mesh("simplified.off", sm, CGAL::parameters::stream_precision(17));

  //alternative simplification; remesh_almost_planar_patches & remove_almost_degenerate_faces 
  //to maintain most of the shape of the original mesh. 
  /*
   bool CGAL::Polygon_mesh_processing::remesh_almost_planar_patches	(	const TriangleMeshIn & 	tm_in,
                                                                                PolygonMeshOut & 	pm_out,
                                                                                std::size_t 	nb_patches,
                                                                                std::size_t 	nb_corners,
                                                                                FacePatchMap 	face_patch_map,
                                                                                VertexCornerMap 	vertex_corner_map,
                                                                                EdgeIsConstrainedMap 	ecm,
                                                                                const NamedParametersIn & 	np_in = parameters::default_values(),
                                                                                const NamedParametersOut & 	np_out = parameters::default_values() 
                                                                        )
  


  */
  //remake surface mesh for testing so we don't use a pre-simplified version
  Tr tr2;
  C2t3 c2_t3(tr2);

  CGAL::make_surface_mesh(c2_t3, surface, criteria, CGAL::Non_manifold_tag());
  Surface_mesh sm2;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2_t3,sm2);
  
  //need typical edge length for isotropic remeshing
  float meanLength = meanEdgeLength(sm2);
  std::cout << "Mean edge length is " << meanLength << std::endl;
  
  Face_range sm2Faces = sm2.faces();
  PMP::isotropic_remeshing(sm2Faces, meanLength, sm2);
  CGAL::IO::write_polygon_mesh("heightmesh_isotropic_remesh.off", sm2, CGAL::parameters::stream_precision(17));
 
  std::cout << "# points after iso remesh " << tr2.number_of_vertices() << "\n";

}
