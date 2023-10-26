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
//types and names
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
typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor  vertex_descriptor;
//typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor    face_descriptor; //maybe necessary but it will be easier to just work with indices
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;


std::vector<Point_3> getVertexPositions(Triangle_mesh mesh, Face_index fIndex) {
   
  Triangle_mesh::Halfedge_index hf = mesh.halfedge(fIndex);
  std::vector<Point_3> vertices; 

  for(Triangle_mesh::Halfedge_index hi : halfedges_around_face(hf, mesh))
  {
    Triangle_mesh::Vertex_index vi = target(hi, mesh);
    vertices.push_back( mesh.point(vi)); //working with xyz points rather than indices --
                                         //don't need to alter base mesh, so points fine
  }
  
  return vertices;
}

std::vector<Vertex_index> getVertexIndices(Triangle_mesh mesh, Face_index fIndex) {
   
  Halfedge_index hf = mesh.halfedge(fIndex);
  std::vector<Vertex_index> vertices; 

  for(Halfedge_index hi : halfedges_around_face(hf, mesh))
  {
    Vertex_index vi = target(hi, mesh);
    vertices.push_back(vi);
  }
  
  return vertices;
}


Face_index getTargetFace(std::vector<Vertex_index> intersected, Point_3 pos, Vector_3 toIntersection, Face_index source_face, Triangle_mesh mesh) {
 //Find the face we're moving into by evaluating the other face attached to the edge made of vertex indices "intersected."

 bool throughVertex = if(intersected.size() = 1);
 if (throughVertex) {
   double moveEpsilon = 1.05; //using a tiny movement in the direction of the intersection vector to determine which face we're moving into
   return PMP::locate(pos+moveEpsilon*toIntersection,mesh).first;//this is really, genuinely, just an approximation so i can debug the rest. 
                                                                 //But it should identify things just fine most of the time.
 }
 
 Halfedge_index intersected_edge = mesh.halfedge(intersected[0],intersected[1]);
 Face_index provisional_target_face = mesh.face(intersected_edge);
 if (provisional_target_face == source_face) provisional_target_face = mesh.face(mesh.opposite(intersected_edge));
 return provisional_target_face;
 
}



Point_3 shift(Triangle_mesh mesh, const Point_3 pos, const Vector_3 move) {
  double travelLength = vectorMagnitude(move);
   
  //we will eventually draw all vertex/edge segments of current face and store in lists; 
  Face_location sourceLocation = PMP::locate(pos, mesh);
  Point_3 source_point = PMP::construct_point(sourceLocation, mesh); //xyz point representing current source
  Vector_3 current_move = move; 
  Face_index currentSourceFace = sourceLocation.first;
  //initializations
  std::vector<Point_3> vertexList;
  std::vector<Segment_3> edgesList; 
  Face_index currentTargetFace;
  Vector_3 currentTargetFaceNormal;
  double lengthToSharedElement;
  std::vector<Face_index> faceIndexList; //store face indices here so we know where to look later

  //useful items for loop w/definition
  bool intersection = true; // true until we have checked all the edges/vertex and verified there's no intersection
  std::size_t counter = 0;
  //find_intersection_baryroutine(Point_3 source, Point_3 target,  std::vector<Point_3> faceVertices) finds the intersection point on an edge
  
  const std::string intersections_filename = "intersections.txt";
  const std::string vertices_filename = "faces_visited.txt";
  std::ofstream intersections_file(intersections_filename);
  std::ofstream vertices_file(vertices_filename);
  
  while(intersection){
    counter += 1;
    std::cout << "" << std::endl;
    std::cout << "ITERATION: " << counter << std::endl;
    std::cout << "current source face is " << currentSourceFace << std::endl; 
    vertexList = getVertexPositions(mesh,currentSourceFace);
    edgesList = createEdgeSegments(vertexList);
    //may need to add a check to see if we're going through a vertex later when finding target face 
    std::cout << "current vertices are: " << std::endl; 
    for (Point_3 vert: vertexList) std::cout << vert << std::endl;  

    //write to files for visualization
    if (intersections_file.is_open()){
      intersections_file << "{" << source_point.x() << ", " << source_point.y() << ", " << source_point.z() << "}";
      intersections_file << "\n";
    }
    if (vertices_file.is_open()) {
      for (Point_3 vert: vertexList) vertices_file << "{" << vert.x() << ", " << vert.y() << ", " << vert.z() << "}" << "\n";
      vertices_file << "\n";
    }

    Point_3 intersection_point = find_intersection_baryroutine(source_point, source_point+current_move, vertexList); 

    if (intersection_point == Point_3(1000,1000,1000)) {
      intersection = false;
      std::cout << "no intersection found after " << counter << " iterations." << std::endl; 
      break; 
    }

    std::cout << "intersection point at " << intersection_point << std::endl;
    Vector_3 vector_to_intersection = Vector_3(source_point, intersection_point);
    double lengthToSharedElement = vectorMagnitude(vector_to_intersection); //how far we've traveled
    current_move = reduceMove(current_move, lengthToSharedElement);         //decrease move size by length to intersected vertex/edge --
                                                                            //effectively the step where we "walk" to that intersection


    currentTargetFace = getTargetFace(source_point, vector_to_intersection, currentSourceFace, mesh); //face we're about to walk into;

    source_point = intersection_point;//update source to be the most recent intersection point, as we have finished walking there

    currentTargetFaceNormal = PMP::compute_face_normal(currentTargetFace,mesh);
    std::cout << "bending the path" << std::endl;

    current_move = projectMoveDown(source_point, currentTargetFaceNormal, normalizer(current_move), vectorMagnitude(current_move)); //bend the path we are about to walk into the plane of the current face
    std::cout<< "Current move is " << current_move << std::endl;

    currentSourceFace = currentTargetFace;
  }
  //source_point+move is the location in the original face if there were no intersections, and it will 
  //be the location in the unfolded mesh if there were intersections (from an edge intersection to a spot
  //within a face)  
  Point_3 fTarget = source_point+current_move;
  intersections_file << "{" << fTarget.x() << ", " << fTarget.y() << ", " << fTarget.z() << "}";
  intersections_file.close();
  vertices_file.close();
  
  return fTarget; //might need some locating/more closely tying this to the mesh, but this should in principle be correct 

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

 Vector_3 forceDisplacement = -20*Vector_3(-0.037457, 0.0330185, 0.177704); //reversing direction because it looks more promising for tests when plotting everything in mathematica 
 Point_3 pointToMove = Point_3(3.51033, 1.9177, 0);


 Point_3 newPos = shift(tmesh, pointToMove, forceDisplacement);
 std::cout << "The post-shift XYZ location is " << newPos << std::endl;
 return 0;
}