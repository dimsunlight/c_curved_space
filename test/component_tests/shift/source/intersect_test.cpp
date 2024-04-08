
//Author: Toler H. Webb
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

 std::vector<Face_index> fs;
 for(Face_index fd: tmesh.faces()){
   fs.push_back(fd);
 }

 Face_index problemFace = fs[284];
 std::cout << "checking against face " << problemFace << std::endl; // ctrl+r redoes undone changes
 std::vector<Point_3> problemVerts = getVertexPositions(tmesh, problemFace); 
 std::vector<Vertex_index> problemIndices = getVertexIndices(tmesh, problemFace); 
 std::cout << "with vertices "; 
 for (Point_3 v: problemVerts) std::cout << v << std::endl; 

 /* barycentric test info: 
  * source point was 0.0283666, 0.971633, 1.45668e-16
  * target point was 0.0152188, 1.03581, -0.0510287
  * intersection point was 0.971633, 7.28357e-18, 0.0283666
  */
 
 std::array<double,3> spointBary = {0.0283666, 0.971633, 1.45668e-16};
 std::array<double,3> tpointBary = {0.0152188, 1.03581, -0.0510287};
 Point_3 spointR3 = PMP::construct_point(std::make_pair(problemFace,spointBary),tmesh);
 Point_3 tpointR3 = PMP::construct_point(std::make_pair(problemFace,tpointBary),tmesh);
 
 std::pair<Point_3,std::vector<Vertex_index>> intersectionInfo = find_intersection(tmesh, problemFace, spointR3, tpointR3,  problemIndices);
 
 std::cout << "intersection point " << intersectionInfo.first << std::endl;
 std::cout << "intersected elements: ";
 for (Vertex_index vi: intersectionInfo.second) std::cout << vi << " "; 


 return 0;
}
