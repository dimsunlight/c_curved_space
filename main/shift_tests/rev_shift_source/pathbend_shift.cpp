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
typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor  vertex_descriptor;
typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor    face_descriptor;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;

//utility functions
auto vectorMagnitude(Vector_3 v) {
  auto const slen = v.x()*v.x() + v.y()*v.y()+v.z()*v.z();
  auto d = CGAL::approximate_sqrt(slen);
  return d; 
}

auto vectorMagnitude(Vector_2 v) {
  auto const slen = v.x()*v.x() + v.y()*v.y();
  auto d = CGAL::approximate_sqrt(slen);
  return d;
}

auto normalizer(Vector_3 v)
{
  auto const d = vectorMagnitude(v);
  return v / d;
}

auto normalizer(Vector_2 v) 
{
  auto const d = vectorMagnitude(v);
  return v/d;
}

int findIndex(Point_3 loc,std::vector<Point_3> vec) {
  auto it = find(vec.begin(), vec.end(), loc);
  int index;
  if (it != vec.end()) {
    index = it - vec.begin();
  }
  else {
    std::cout << "could not find index for entry (" << loc << ")" << std::endl;
    index = -1; 
  }
  return index;
}  


//function definitions
std::vector<Point_3> getVertexPositions(Triangle_mesh mesh, face_descriptor fIndex) {
   
  Triangle_mesh::Halfedge_index hf = mesh.halfedge(fIndex);
  std::vector<Point_3> vertices; //I think we can say 3 because we know it's a triangle mesh

  for(Triangle_mesh::Halfedge_index hi : halfedges_around_face(hf, mesh))
  {
    Triangle_mesh::Vertex_index vi = target(hi, mesh);
    vertices.push_back( mesh.point(vi)); //working with xyz points rather than indices --
                                         //don't need to alter base mesh, so points fine
  }
  
  return vertices;
}

auto getSharedElements(std::vector<Point_3> vs1, std::vector<Point_3> vs2) {
  //common element finder for shared edge; double for loop is slow, but
  //that's fine given guaranteed small size. We could use set intersection finder, 
  //but that needs a sorted vector and I can't sort Point_3 objects!
  std::vector<Point_3> sharedElements; 
  for (Point_3 i: vs1) {
    for (Point_3 j: vs2) {
      if (i == j) {
        sharedElements.push_back(i);
      }
    }
  }
  if (sharedElements.size() == 1) std::cout << "Only 1 vertex shared." << std::endl;
  if (sharedElements.size() == 0) std::cout << "No shared vertices." << std::endl;
  return sharedElements;
}

auto getUnsharedElements(std::vector<Point_3> vertices, std::vector<Point_3> shared) {
  //find vertices not in "shared."  
  //Shared is the vector of vertex locations
 
  std::vector<Point_3> oddElements;
  
  for (Point_3 i: vertices) {
    if (std::find(shared.begin(), shared.end(), i) != shared.end()) {}
    else {
      oddElements.push_back(i);
    }
  }

  return oddElements;
}

double angleBetween(face_descriptor f1, face_descriptor f2, Triangle_mesh mesh) {
  //generate normals -- compute_face_normal gives unit normals by default
  Vector_3 oldNormal = PMP::compute_face_normal(f1,mesh), newNormal = PMP::compute_face_normal(f2,mesh);
  //normals compute fine (verified via print) 

  FT oldDotNew = CGAL::scalar_product(oldNormal,newNormal);
  double angle = acos(oldDotNew);
  //angle within expectations unit test 
  if (angle > 1) {
    std::cout << "Angle between normals over 1!" << std::endl;
  }

  return angle;
}

auto rotateAboutAxis(std::vector<Point_3> targets, std::vector<Point_3> axis, double rotAngle) {
  //rotation about an axis named axis; handles rotation of multiple targets at once
  
  Vector_3 axisVector = normalizer(Vector_3(axis[0],axis[1]));
  std::vector<Point_3> shiftedTargets; 

  for (Point_3 target: targets) {
    shiftedTargets.push_back(target - Vector_3({0,0,0}, axis[0]));
  }
  std::cout << "targets shifted for rotation to " << std::endl;
  for (Point_3 shiftedTarget: shiftedTargets) std::cout << shiftedTarget << std::endl;

  double u1 = axisVector.x(), u2 = axisVector.y(), u3=axisVector.z();
  
  double c = cos(rotAngle), s = sin(rotAngle);
  double C = 1-c;
  //unfortunately, I believe I have to multiply element wise in C++
  double Identity[3][3] = {{c,0.0,0.0},{0.0,c,0.0},{0.0,0.0,c}};
  double crossMatrix[3][3] = {{0,-s*u3,s*u2},{s*u3,0,-s*u1},{-s*u2,s*u1,0}};
  double tensorProduct[3][3] = {{C*u1*u1,C*u1*u2,C*u1*u3},{C*u2*u1,C*u2*u2,C*u2*u3},{C*u3*u1,C*u3*u2,C*u3*u3}};
  
  double rMatrix[3][3];
  //this double loop to define the rotation matrix is probably a place where it'd be easy to speed things up very slightly
  //with more c++ knowledge
  for (std::size_t i = 0; i < 3; i++) {
    for (std::size_t j = 0; j < 3; j++) {
	    rMatrix[i][j] = Identity[i][j] + crossMatrix[i][j] + tensorProduct[i][j];
    }
  }
  
  std::vector<Point_3> rotatedTargets;
  double rotatedStorage[3];
  std::vector<double> targetStorage;
  
  //rotate and shift back all targets
  for (Point_3 shiftedTarget: shiftedTargets) {
    targetStorage = {shiftedTarget.x(),shiftedTarget.y(),shiftedTarget.z()};
  
    for (std::size_t i = 0; i < 3; i++) {
      double product = 0;
      for (std::size_t j = 0; j < 3; j++) {
        product += rMatrix[i][j]*targetStorage[j];
      }
      rotatedStorage[i] = product; 
    }
    Point_3 rotatedTarget = Point_3(rotatedStorage[0],rotatedStorage[1],rotatedStorage[2]);
    
    //undo shift at the end
    rotatedTargets.push_back(rotatedTarget + Vector_3({0,0,0}, axis[0]));
  }

  return rotatedTargets;
}

auto createTemporaryMesh(std::vector<Point_3> vertexTrio) {
  Triangle_mesh tempMesh;
  vertex_descriptor u = tempMesh.add_vertex(vertexTrio[0]);
  vertex_descriptor v = tempMesh.add_vertex(vertexTrio[1]);
  vertex_descriptor w = tempMesh.add_vertex(vertexTrio[2]);
  face_descriptor tFace = tempMesh.add_face(u,v,w);
  std::cout << "temp face normal = " << PMP::compute_face_normal(tFace,tempMesh) << std::endl;
  return tempMesh;
}

Point_3 getNewXYZLocation(Point_3 flatLocation, Triangle_mesh originalMesh, Triangle_mesh tempMesh, face_descriptor f2) {
  Face_location tempLocation= PMP::locate(flatLocation,tempMesh);//get bary location in temporary mesh
  Barycentric_coordinates b_weights = tempLocation.second;//export bary location as set of bary weights
  Face_location frankenLocation = std::make_pair(f2, b_weights); //create location object from original face index + bary weights
  Point_3 newXYZLocation = PMP::construct_point(frankenLocation, originalMesh); //use new location object to find location in original mesh
  return newXYZLocation;
}


std::vector<Segment_3> createEdgeSegments(std::vector<Point_3> vs) {
  //create Segment_3 objects to match all the edges of a face. 
  std::vector<Segment_3> edges;
  std::size_t numVertices = vs.size();//size is 1 more than the highest index

  for (std::size_t i; i < numVertices; i++) {
    std::size_t j = i+1; 
    while(j < numVertices) {
      edges.push_back(Segment_3(vs[i],vs[j]));
      j += 1;
    }
  }
  //std::cout << " printing edges in edge list" << std::endl;
  //for (Segment_3 edge: edges) std::cout << edge << std::endl;


  return edges;
}



face_descriptor getTargetFace(Point_3 pos, Vector_3 toIntersection, face_descriptor currentSourceFace, Triangle_mesh mesh) {
 //goal is to find the face we're moving into. this is very easy if we share an edge, and harder if we only share a vertex. 
 //easy if we share an edge because we can just look for faces which share those vertices;
 //hard if we only share a vertex because the faces that share a vertex are not limited to the source and target.
 
 //big search function will be too slow... maybe?
 double moveEpsilon = 1.05; //using a tiny movement in the direction of the intersection vector to determine which face we're moving into

 return PMP::locate(pos+moveEpsilon*toIntersection,mesh).first;//this is really, genuinely, just an approximation so i can debug the rest. 
                                                               //But it should identify things just fine most of the time. 
}

Vector_3 reduceMove(Vector_3 moveVector, double reduceLengthBy) {
  double originalLength = vectorMagnitude(moveVector);
  Vector_3 normalizedMove = normalizer(moveVector); 
  Vector_3 reducedMove = (originalLength-reduceLengthBy)*normalizedMove;

  return reducedMove;
}

Vector_3 projectMoveDown(Point_3 source, Vector_3 targetFaceNormal, Vector_3 moveDirection, double remainingDistance) {
  //takes as input a source point, face normal, movement direction vector, and distance left to travel,
  //and returns a new movement vector from the source point along the surface of the current face. 
  double overlap = moveDirection * targetFaceNormal; //this is inner product in CGAL if both are Vector_3 objects
  Vector_3 inFaceVector = remainingDistance*normalizer(moveDirection - overlap*targetFaceNormal);
  return Vector_3(source, source+inFaceVector); 
}

double triangle_area(Point_3 v1, Point_3 v2, Point_3 v3) {
  Vector_3 side1 = Vector_3(v1, v2);
  Vector_3 side2 = Vector_3(v1, v3);
  return vectorMagnitude(CGAL::cross_product(side1,side2))/2;
}

Point_3 project_to_face(std::vector<Point_3> vertices, Point_3 query) {
  Point_3 projectedQuery;
  Vector_3 side1, side2, normal, toQuery;
  side1 = Vector_3(vertices[0],vertices[1]);
  side2 = Vector_3(vertices[0],vertices[2]);
  normal = CGAL::cross_product(side1, side2); //don't need to worry about direction due to dot
  toQuery = Vector_3(vertices[0], query);
  projectedQuery = vertices[0] + (toQuery-CGAL::scalar_product(toQuery,normal)*normal);
  return projectedQuery;
}


Point_3 find_intersection_baryroutine(Point_3 source, Point_3 target,  std::vector<Point_3> faceVertices) {
  //use the barycentric coordinates of a point in order to determine the R3 coordinates of its intersection with an edge of the current face, if such an intersection exists.
  //takes only the r3 positions of the source and prospective target along with the vertices of the current face; intersection test exists independent of the mesh.  
  Point_3 sourceDown = project_to_face(faceVertices, source);
  Point_3 query = project_to_face(faceVertices, target);
  std::array<double, 3> source_bary = PMP::barycentric_coordinates(faceVertices[0],faceVertices[1],faceVertices[2], sourceDown);
  std::array<double, 3> query_bary = PMP::barycentric_coordinates(faceVertices[0],faceVertices[1],faceVertices[2], query);
  Vector_3 source_bary_vector = Vector_3(source_bary[0],source_bary[1],source_bary[2]); //define as vectors so we can manipulate them via addition/subtraction/etc
  Vector_3 query_bary_vector  = Vector_3(query_bary[0], query_bary[1] , query_bary[2]);
  Vector_3 displacement = query_bary_vector - source_bary_vector;
  //defining bary components as constants for easy manipulation
  double const b11 = source_bary_vector[0];
  double const b12 = source_bary_vector[1];
  double const b13 = source_bary_vector[2];
  //each entry of intersection_values is the first point along the coordinate component of the line in that barycentric coordinate
  //where the barycentric coordinate reaches 0; when a barycentric coordinate reaches 0, we've intersected an edge. The first one
  //(the smallest intersection_value) is the first time we hit an edge, and should therefore be the intersection point.The others
  //will indicate how far we'd need to travel along the displacement vector before another coordinate became 0 (and we'd intersect a "phantom" edge)
  std::vector<double> intersection_values = {-b11/displacement[0], -b12/displacement[1], -b13/displacement[2]};
  //Entries are positive unless we'd never make a barycentric weight 0 by traveling along the displacement; we discard the negative
  //results which correspond to that.
  double toIntersect = intersection_values[0];
  double tol = 0.0001;
  //the first value is the only one not checked against tol in the main minimum-value-finding routine, so. 
  //we check it now. 
  if (toIntersect < tol) toIntersect = 2;

  //min function with a max of 1. if we don't find something less than 1, no intersection.  
  std::cout << " calling intersection finder. " << std::endl;
  std::cout << "intersection values are: " << std::endl;

  for(double val: intersection_values) {
    std::cout << val << std::endl;
    if (val < toIntersect and val > tol) toIntersect = val;

  }

  std::cout << "toIntersect is " << toIntersect << std::endl;
  if (toIntersect < 1 and toIntersect > tol) {
    std::cout << "found intersection. " << std::endl;
    Point_3 min_intersection;
    min_intersection = Point_3(b11,b12,b13);
    min_intersection = min_intersection+toIntersect*displacement;
    std::cout << "barycentric intersection point " << min_intersection << std::endl;
    Vector_3 xyz_intersection = min_intersection[0]*Vector_3(faceVertices[0].x(),faceVertices[0].y(),faceVertices[0].z()) + min_intersection[1]*Vector_3(faceVertices[1].x(),faceVertices[1].y(),faceVertices[1].z())+ min_intersection[2]*Vector_3(faceVertices[2].x(),faceVertices[2].y(),faceVertices[2].z()); //construct point from definition of barycentric to xyz conversion
    std::cout << "xyz intersection point (in f_i) " << xyz_intersection << std::endl;
    Point_3 intersection_Point_3 = Point_3(0,0,0) + xyz_intersection;
    return intersection_Point_3; // this version returns the barycentric intersection point
  }
  else {
    std::cout << "No intersection. Returning filler point.";
    return Point_3(1000,1000,1000);
  }
  return Point_3(10000,1000,1000); //default return
}

Point_3 shift(Triangle_mesh mesh, const Point_3 pos, const Vector_3 move) {
  double travelLength = vectorMagnitude(move);
   
  //we will eventually draw all vertex/edge segments of current face and store in lists; 
  Face_location sourceLocation = PMP::locate(pos, mesh);
  Point_3 source_point = PMP::construct_point(sourceLocation, mesh); //xyz point representing current source
  Vector_3 current_move = move; 
  face_descriptor currentSourceFace = sourceLocation.first;
  std::cout << "source point after locate at " << source_point << std::endl;
  //initializations
  std::vector<Point_3> vertexList;
  std::vector<Segment_3> edgesList; 
  face_descriptor currentTargetFace;
  Vector_3 currentTargetFaceNormal;
  double lengthToSharedElement;
  std::vector<face_descriptor> faceIndexList; //store face indices here so we know where to look later

  //useful items for loop w/definition
  bool intersection = true; // true until we have checked all the edges/vertex and verified there's no intersection
  std::size_t counter = 0;
  //find_intersection_baryroutine(Point_3 source, Point_3 target,  std::vector<Point_3> faceVertices) finds the intersection point on an edge
  
  const std::string intersections_filename = "intersections.txt";
  const std::string vertices_filename = "faces_visited.txt";
  std::ofstream intersections_file(intersections_filename);
  std::ofstream vertices_file(vertices_filename);
  
  while(intersection){
    std::cout << "" <<std::endl;
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


