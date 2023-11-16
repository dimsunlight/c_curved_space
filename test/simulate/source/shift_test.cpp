
//Author: Toler H. Webb
//Description: code which takes as input a position on a mesh and a 
//displacement and returns as output a new position on the mesh. 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <math.h> //including this to do simple cos and sine, but maybe imprecise compared to cgal-originating routines
#include "utils.h" 
#include "shift.h"
#include <chrono>
#include <ratio> 

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::FT                                                           FT;
typedef K::Point_2                                                      Point_2;
typedef K::Ray_2                                                        Ray_2;
typedef K::Point_3                                                      Point_3;
typedef K::Vector_3                                                     Vector_3;
typedef K::Ray_3                                                        Ray_3;
typedef K::Segment_3                                                    Segment_3;
typedef K::Intersect_3                                                  Intersect_3;
typedef K::Triangle_3                                                   Triangle_3;
typedef K::Segment_2                                                    Segment_2;
typedef K::Point_2                                                      Point_2;
typedef K::Vector_2                                                     Vector_2;
typedef CGAL::Surface_mesh<Point_3>                                     Triangle_mesh;
typedef Triangle_mesh::Face_index                                       Face_index;
typedef Triangle_mesh::Halfedge_index                                   Halfedge_index;
typedef Triangle_mesh::Vertex_index                                     Vertex_index;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;

float meanTriangleArea(Triangle_mesh mesh) {
  float mean = 0;
  int count = 0;
  auto facelist = mesh.faces();
  for (Face_index face: facelist) {
    ++count;
    std::vector<Point_3> vs = getVertexPositions(mesh,face);
    mean += triangle_area(vs[0],vs[1],vs[2]);
  }
  return mean/count;
}

int main(int argc, char* argv[]) {

 //get input mesh from command line argument
 const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("sims_project/torusrb20.off");
 Triangle_mesh tmesh;
 if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
   !CGAL::is_triangle_mesh(tmesh))
 {
   std::cerr << "Invalid input file." << std::endl;
   return EXIT_FAILURE;
 }
 
 float sqrtMTA = sqrt(meanTriangleArea(tmesh)); 
 //using my torus source shift... i should probably go into shift and make sure
 //this is forced tangent to the face in question
 Vector_3 forceDisplacement = normalizer(Vector_3(-0.037457, 0.0330185, 0.177704)); //reversing direction because it looks more promising for tests when plotting everything in mathematica

 //micro-projection routine
 Point_3 dummyPt = Point_3(3.51033, 1.9177, 0); //intended source point
 Face_location dummyLoc = PMP::locate(dummyPt, tmesh); //finding where it would be on the mesh
 dummyPt = PMP::construct_point(dummyLoc, tmesh); //getting the closest point to the intended point that's on the mesh
 std::vector<Point_3> sFaceVerts = getVertexPositions(tmesh,dummyLoc.first); 
 Point_3 newTarget = project_to_face(sFaceVerts,dummyPt+forceDisplacement);
 Vector_3 dummyDisplacement = normalizer(Vector_3(dummyPt,newTarget)); 

 //now loop over a few sizes relative to the average length scale of the mesh and see how long it takes 
 std::string times_filename = filename.substr(0,9) + "shift.csv"; 
 std::ofstream times_file(times_filename);
 times_file << filename << ", Mean Triangle Area " << meanTriangleArea(tmesh) << "\n"; 
 

 std::chrono::duration<long, std::milli> s_time;

 for (int i = 0; i < 20; i++) {
   std::cout << "iteration " << i << " shift length " << i*sqrtMTA <<  std::endl;
   auto start = std::chrono::high_resolution_clock::now();
   Point_3 np = shift(tmesh, dummyPt, i*sqrtMTA*dummyDisplacement);
   auto end = std::chrono::high_resolution_clock::now();
   s_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
   times_file << i*sqrtMTA << ", " << s_time.count() << "\n"; 
 }

 return 0;
}
