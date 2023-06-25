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
//types and names
typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::FT                                                           FT;
typedef K::Point_2                                                      Point_2;
typedef K::Ray_2                                                        Ray_2;
typedef K::Point_3                                                      Point_3;
typedef K::Vector_3                                                     Vector_3;
typedef K::Ray_3                                                        Ray_3;
typedef CGAL::Surface_mesh<Point_3>                                     Triangle_mesh;
typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor           vertex_descriptor;
typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor             face_descriptor;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                                    Face_location;

//function definitions
auto overEdge(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
  return Point_3(0.0,0.0,0.0);
}

auto shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
  //first: check if we will be in our out of triangle after naive shift
  //then: do the simple shift if we will, unfold shift if we won't 

  Point_3 trialNewPos = pos + move;
  Face_location oldPosLocation = PMP::locate(pos, mesh);
  std::cout << oldPosLocation.first << std::endl;

  Face_location newPosLocation = PMP::locate(trialNewPos,mesh);
  std::cout << newPosLocation.first << std::endl;

  bool withinTriangle=true;

  if (withinTriangle) {
    //could also return bary coordinates; PMP::locate(trialNewPos,mesh)
    return trialNewPos;
  }

  if (!withinTriangle) {
    return overEdge(mesh, pos, move);
  }

  std::cout << "neither in or out of triangle...?" << std::endl;
  return Point_3(0.0,0.0,0.0);
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

 return 0;
}





