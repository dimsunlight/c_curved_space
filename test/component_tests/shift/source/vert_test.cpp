
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
typedef CGAL::Face_around_target_circulator<Triangle_mesh>              Face_circulator;

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
 //	1) vertex & vertex location (to be intersected on purpose) 
 //	2) all the faces around said vertex 
 //	3) create outside a mathematica drawing with vertex, faces, and 
 //	   path information to port back in here! 
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
 return 0;
}
