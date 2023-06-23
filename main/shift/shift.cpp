//Author: Toler H. Webb
//Description: code which takes as input a position on a mesh and a 
//displacement and returns as output a new position on the mesh. 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>


//typedefs:
typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     KernelWithSqrt;
typedef Kernel::Point_3                                                 Point_3;
typedef Kernel::Vector_3                                                Vector_3;
typedef Kernel::Ray_3                                                   Ray_3;
typedef CGAL::Surface_mesh<Point_3>                                     Triangle_mesh;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh>  Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits>                        Surface_mesh_shortest_path;
typedef boost::graph_traits<Triangle_mesh>                              Graph_traits;
typedef Graph_traits::vertex_iterator                                   vertex_iterator;
typedef Graph_traits::face_iterator                                     face_iterator;
typedef typename Surface_mesh_shortest_path::Barycentric_coordinates    Barycentric_coordinates;
typedef typename Surface_mesh_shortest_path::Face_location              Face_location;
typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh>         AABB_face_graph_primitive;
typedef CGAL::AABB_traits<Kernel, AABB_face_graph_primitive>            AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                         AABB_tree;

//namespaces
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;


//function definitions
auto simpleShift(Triangle_mesh mesh,Point_3 pos, Vector_3 move) {
  
  return 0;
}

auto overEdge(Triangle_mesh mesh, point_3 pos, Vector_3 move) {
  
  return 0;
}

auto shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
  //first: check if we will be in our out of triangle after naive shift
  //then: do the simple shift if we will, unfold shift if we won't
  
  const Point_3 query(1.2, 7.4, 0);
  Face_location query_location = PMP::locate(query, tm);

  Point_3 trialNewPos = pos + move;

  
  if (withinTriangle) {
    return simpleShift(mesh, pos, move);
  }

  if (outsideTriangle) {
    return overEdge(mesh, pos, move);
  }

  std::cout << "neither in or out of triangle...?" << std::endl;
  return 0;
}


//main
int main(int argc, char* argv[]) {
 
 //get input mesh from command line argument
 const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("sims_project/torusrb12.off");
 Triangle_mesh tmesh;
 if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
   !CGAL::is_triangle_mesh(tmesh))
 {
   std::cerr << "Invalid input file." << std::endl;
   return EXIT_FAILURE;
 }

 Vector_3 forceDisplacement = Vector_3(-0.0145553, 0.0453073, 0.176461);
 Point_3 pointToMove = Point_3(3.51033, 1.9177, 0);

 Point_3 newPos = shift(tmesh, pointToMove, forceDisplacement);
 std::cout << newPos << std::endl;

 return 0;
}





