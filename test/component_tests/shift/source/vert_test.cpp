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
 //now -- build out mesh information necessary to create a test case. 
 //	1) generate a random bary  point on one of the faces adjoining vertex v.  
 //	2) get the vector from the random point to the vertex location. 
 //	3) test throughvertex from the barypoint to the vertex exactly. 
 //	for each vertex.
 Point_3 vpoint;
 double b1; 
 double b2; 
 double b3; 

 Point_3 rand_source;
 Point_3 rand_end; 

 Vector_3 to_v;
 MTgenerator gen = make_generator(1997); 
 std::vector<Face_index> v_faces; 
 v_faces.reserve(6); 
 Face_index found_target; 
 std::array<double,3> rand_weights;
 Point_3 shift_output;
 double to_dist; 
 printf("starting vertex test loop...\n"); 
 for (Vertex_index v: vs) { 
   vpoint = tmesh.point(v); 
   Face_circulator fbegin(tmesh.halfedge(v),tmesh), done(fbegin); 
   v_faces.clear(); 
   do {
     v_faces.push_back(*fbegin++);
   } while(fbegin != done);
   //defining weights outside of loop/array definition explicitly for later
   b1 = rand_weight(gen); 
   b2 = rand_weight(gen); 
   b3 = rand_weight(gen); 
   rand_weights[0] = b1; 
   rand_weights[1] = b2; 
   rand_weights[2] = b3; 
   Face_location rand_loc = std::make_pair(v_faces[0], rand_weights);
   rand_source = PMP::construct_point(rand_loc, tmesh); 
   to_v = Vector_3(rand_source, vpoint); 
   rand_end = rand_source + 2*to_v;   
 
   double bigShiftMod = 15; // make this number >> 1 to test against very large shifts  

   std::cout << "selectFaceFromVertex for vertex " << v << "..." << std::endl; 
   found_target = selectFaceFromVertex(v, to_v, v_faces[0], tmesh).first;
   std::cout << "For vertex " << v << ":" << std::endl;
   std::cout << "Found target face at: " << found_target << std::endl;
   std::cout << "---STARTING SHIFT---" << std::endl; 
   shift_output = shift(tmesh, rand_source, bigShiftMod*to_v);   
   std::cout << "Shifted to: " << shift_output << std::endl;
   std::cout << "\n"; 

 
   if (v_faces[0] == found_target) {
     //need to see vertex positions for the faces and whatnot 
     std::cout << "Faces: " << std::endl;
     for (Face_index fIndex: v_faces) {
       std::vector<Point_3> vertPos = getVertexPositions(tmesh, fIndex);
       std::cout << fIndex << ": "; 
       for (Point_3 vert: vertPos) {
         std::cout << "{" << vert.x() << ", " << vert.y() << ", " << vert.z() << "}, ";      
       }
       std::cout << std::endl;
     }
     std::cout << "vertex point " << vpoint << std::endl;
     std::cout << "source point " << rand_source << std::endl;   
     std::cout << "Source face: " << v_faces[0] << std::endl; 
     break; 
   }
   if (v == vs[414]) break; 
 }

 return 0;

}
/*
 Vertex_index intersectedVertex = vs[40];
 
 std::cout << "Vertex location: " << tmesh.point(intersectedVertex) << std::endl;

 Face_circulator fbegin(tmesh.halfedge(intersectedVertex),tmesh), done(fbegin); //fbegin? (lol)
 std::vector<Face_index> facesToPrint;
 facesToPrint.reserve(6);
 do {
    facesToPrint.push_back(*fbegin++);
 } while(fbegin != done);
 
 std::vector<Point_3> vertPos; 
 vertPos.reserve(3); 
 
 std::cout << "Faces: " << std::endl;
 for (Face_index fIndex: facesToPrint) {
   vertPos = getVertexPositions(tmesh, fIndex);
   for (Point_3 vert: vertPos) {
     std::cout << "{" << vert.x() << ", " << vert.y() << ", " << vert.z() << "}, ";      
   }
   std::cout << std::endl;
 }

 // having gone into mathematica and figured out a point within a face (src_point),
 // then moved it through the vertex (mid_point), to end up at a target (end_point) 
 // for rotation... 
 Point_3 src_point = Point_3(-2.45268, -0.523961, 0.853549);
 Point_3 mid_point = tmesh.point(intersectedVertex);
 Point_3 end_point = Point_3(-2.30864, -0.461563, 0.784361);
 
 printf("\n");
 std::cout << "Correct faces are: " << facesToPrint[1] << " and " << facesToPrint[4] << std::endl;
 std::cout << "(second and fifth vertex trios from above)" << std::endl;
 
 Face_index newTarget = selectFaceFromVertex(intersectedVertex, Vector_3(mid_point,end_point), facesToPrint[1], tmesh); 

 std::cout << "Routine found the correct face to be: " << newTarget << std::endl;
 */

