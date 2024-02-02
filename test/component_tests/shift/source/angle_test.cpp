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
#include <random> 

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::FT                                                           FT;
typedef K::Point_2                                                      Point_2;
typedef K::Ray_2                                                        Ray_2;
typedef K::Point_3                                                      Point_3;
typedef K::Vector_3                                                     Vector_3;
typedef K::Ray_3                                                        Ray_3;
typedef K::Segment_3                                                    Segment_3;
typedef K::Intersect_3                                                  Intersect;
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
typedef CGAL::Face_around_target_circulator<Triangle_mesh>              Face_circulator;
typedef std::mt19937                                                    MTgenerator; //mersenne twister generator
typedef std::uniform_real_distribution<double>                           distribution;

MTgenerator make_generator(int seed)
{
  //create random number generator using seed parameter. 

  std::random_device                  rand_dev;
  MTgenerator                         generator(rand_dev());
  generator.seed(seed); //seeded for reproducibility -- make sure we use the same sequence every time :)  

  return generator;
}

distribution distribution_of_range(double range_min, double range_max) {
  distribution distro(range_min, range_max);
  return distro;
}

double rand_weight(MTgenerator& gen) {
  //pass by address above to ensure that we sequentially sample
  //from the generator without using the same starting point.
  distribution distro = distribution_of_range(.001, 1); //no points on edges or verts
  double b;
  // reference gen by address so each call to gen advances the
  // rng forward once
  return distro(gen);
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
 
 std::vector<Vertex_index> vs; 
 for(Vertex_index vd: tmesh.vertices()){
   vs.push_back(vd);
 } 

 Vector_3 to_v;
 std::vector<Face_index> v_faces; 
 v_faces.reserve(6); 
 Face_index found_target; 
 Point_3 shift_output;

 Vertex_index intersectedVertex = vs[414];
 Point_3 vPoint = tmesh.point(intersectedVertex); 

 std::cout << "Vertex location: " << vPoint << std::endl;
 std::cout << "vertex id: " << intersectedVertex << std::endl;
 Face_circulator fbegin(tmesh.halfedge(intersectedVertex),tmesh), done(fbegin); 
 std::vector<Face_index> facesToPrint;
 facesToPrint.reserve(6);
 do {
    facesToPrint.push_back(*fbegin++);
 } while(fbegin != done);
 
 std::vector<Point_3> vertPos; 
 vertPos.reserve(3); 
 
 std::array<double,3> arbibary = {.3,.4,.3}; 

 Face_location testPoint = std::make_pair(facesToPrint[0], arbibary); 
 Point_3 inFace = PMP::construct_point(testPoint, tmesh); 
 to_v = Vector_3(inFace, vPoint); //placeholder to check basic functionality 
 Face_index sface = facesToPrint[0]; 
 
 std::cout << "later, vertex id: " << intersectedVertex << std::endl;
 throughVertexByAngle(intersectedVertex, to_v, sface, tmesh);
 
 return 0;
}
