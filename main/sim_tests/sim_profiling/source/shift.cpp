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
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;

Face_index getTargetFace(std::vector<Vertex_index> intersected, Point_3 pos, Vector_3 toIntersection, Face_index source_face, Triangle_mesh mesh) {
 //Find the face we're moving into by evaluating the other face attached to the edge made of vertex indices "intersected."

 bool throughVertex = intersected.size() == 1;
 if (throughVertex) {
   std::cout << "Only one vertex intersected. Calling placeholder routine." << std::endl;
   double moveEpsilon = 1.05; //using a tiny movement in the direction of the intersection vector to determine which face we're moving into
   return PMP::locate(pos+moveEpsilon*toIntersection,mesh).first;//this is really, genuinely, just an approximation so i can debug the rest. 
                                                                 //But it should identify things just fine most of the time.
 }
 
 Halfedge_index intersected_edge = mesh.halfedge(intersected[0],intersected[1]);
 Face_index provisional_target_face = mesh.face(intersected_edge);
 if (provisional_target_face == source_face) provisional_target_face = mesh.face(mesh.opposite(intersected_edge));
 return provisional_target_face;
}

std::pair<Point_3,std::vector<Vertex_index>> find_intersection(Triangle_mesh mesh, Face_index sourceFace, Point_3 source, Point_3 target,  std::vector<Vertex_index> vertexIndices) {
  //use the barycentric coordinates of a point in order to determine the R3 coordinates of its intersection with an edge of the current face, if such an intersection exists.
  //takes only the r3 positions of the source and prospective target along with the vertices of the current face; intersection test exists independent of the mesh.
  
  //need to update this to find the intersected edge. (as a vertex pair) 
  std::vector<Point_3> faceVertices;
  for (Vertex_index vi: vertexIndices) faceVertices.push_back(mesh.point(vi));
   
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
  double tol = 0.0001; //checking for equivalence to zero within reasonable error

  //the first value is the only one not checked against tol in the main minimum-value-finding routine, so. 
  //we check it now. 
  if (toIntersect < tol) toIntersect = 2;//just needs to be some # greater than 1. 

  //min function with a max of 1. if we don't find something less than 1, no intersection.  
  for(double val: intersection_values) {
    if (val < toIntersect and val > tol) toIntersect = val;
  }
    
  std::vector<Vertex_index> fillerVector = {}; //for null return later
  
  if (toIntersect < 1 and toIntersect > tol) {
    Point_3 min_intersection;
    min_intersection = Point_3(b11,b12,b13);
    min_intersection = min_intersection+toIntersect*displacement;
    std::array<double,3> newBaryCoords = {min_intersection.x(),min_intersection.y(),min_intersection.z()};
    
    //off vertices tell us what vertices are *not* part of the intersected edge
    int vertexIndicator = 0;
    std::vector<int> offVertices = {}; //indices of vertices with bary coordinate zero at the intersection
    //would use more economical "erase" for the below, but removing entries will mess up the indices if we have to erase twice, 
    //and we don't know a priori which index is being removed/in what order. 
    for (double b: newBaryCoords) {
      if (b < tol) {
        offVertices.push_back(vertexIndicator);
      }
      vertexIndicator+=1;
    }

    std::vector<Vertex_index> intersected = {};
    bool flag = false;
    for (int i; i < vertexIndices.size(); i++) {

      for (int off: offVertices) {
        if (i==off) flag = true;
      }
      if (flag == true) {
        flag = false;
	continue;
      }

      intersected.push_back(vertexIndices[i]);
    }

    Face_location newPosition = std::make_pair(sourceFace,newBaryCoords);
    Point_3 xyz_intersection = PMP::construct_point(newPosition,mesh);
    /*
    //manual reconstruction  
    Vector_3 xyz_intersection = min_intersection[0]*Vector_3(faceVertices[0].x(),faceVertices[0].y(),faceVertices[0].z()) + 
	                        min_intersection[1]*Vector_3(faceVertices[1].x(),faceVertices[1].y(),faceVertices[1].z()) +
			       	min_intersection[2]*Vector_3(faceVertices[2].x(),faceVertices[2].y(),faceVertices[2].z()); 
    */

    return std::make_pair(xyz_intersection,intersected); // this version returns the barycentric intersection point
  }
  //if we don't have an intersection -- FIX CONDITIONALS LATER IN THE CODE
  else {
    std::cout << "No intersection. Returning filler point.";
    return std::make_pair(Point_3(1000,1000,1000),fillerVector);
  }
  return std::make_pair(Point_3(10000,1000,1000),fillerVector); //default return
}

Face_location rotateIntoNewFace(Triangle_mesh mesh, Face_index sface,
	                         	  Face_index tface, Point_3 source, Point_3 target) {
  //returns the barycentric coordinates of the shift vector currently in sface after
  //rotation to the tangent plane of tface.

  Vector_3 sNormal = PMP::compute_face_normal(sface,mesh);
  Vector_3 tNormal = PMP::compute_face_normal(tface,mesh);
  double angle = angleBetween(sface,tface, mesh); 
  Vector_3 axisVector = normalizer(CGAL::cross_product(sNormal,tNormal));
  std::vector<Point_3> axis = {source, source+axisVector};
  Point_3 rotatedT = rotateAboutAxis(target, axis, angle); 
  std::vector<Point_3> tfaceVertices = getVertexPositions(mesh, tface);
  std::array<double, 3> rotated_barycentric = PMP::barycentric_coordinates(tfaceVertices[0],tfaceVertices[1],
		                                                           tfaceVertices[2],rotatedT);
  Face_location rotatedPosition = std::make_pair(tface,rotated_barycentric);
  return rotatedPosition; 
}

Point_3 shift(Triangle_mesh mesh, const Point_3 pos, const Vector_3 move) {
  double travelLength = vectorMagnitude(move);

  //we will eventually draw all vertex/edge segments of current face and store in lists; 
  Face_location sourceLocation = PMP::locate(pos, mesh);
  Point_3 source_point = PMP::construct_point(sourceLocation, mesh); //xyz point representing current source
  Vector_3 current_move = move;
  Face_index currentSourceFace = sourceLocation.first;
  //initializations
  std::vector<Vertex_index> vertexList;
  std::vector<Point_3> targetVertices;
  std::vector<Point_3> sharedEdge;
  std::vector<Point_3> forRotation;
  Point_3 target;
  Point_3 rotatedTarget;
  std::vector<Segment_3> edgesList;
  Face_index currentTargetFace;
  Vector_3 currentTargetFaceNormal;
  double lengthToSharedElement;
  double rotationAngle;
  double overlap;
  std::vector<Face_index> faceIndexList; //store face indices here so we know where to look later

  //useful items for loop w/definition
  bool intersection = true; // true until we have checked all the edges/vertex and verified there's no intersection
  std::size_t counter = 0;
  //find_intersection_baryroutine(Point_3 source, Point_3 target,  std::vector<Point_3> faceVertices) finds the intersection point on an edge

  const std::string intersections_filename = "r_intersections.txt";
  const std::string vertices_filename = "r_faces_visited.txt";
  const std::string cmove_filename = "current_move_list.txt"; 
  std::ofstream intersections_file(intersections_filename);
  std::ofstream vertices_file(vertices_filename);
  std::ofstream cmove_file(cmove_filename); 

  //assuming we've fed in a vector tangent to the source face, we can use
  //its level of normal overlap as the baseline of normal overlap for rotations
  Vector_3 currentSourceNormal = PMP::compute_face_normal(currentSourceFace,mesh);
  double baseAgreement = currentSourceNormal*move;

  while(intersection){
    counter += 1;
    std::cout << "" << std::endl;
    std::cout << "ITERATION: " << counter << std::endl;
    //vertexList = getVertexPositions(mesh,currentSourceFace);
    vertexList = getVertexIndices(mesh,currentSourceFace);

    if (intersections_file.is_open()){
      intersections_file << "{" << source_point.x() << ", " << source_point.y() << ", " << source_point.z() << "}";
      intersections_file << "\n";
    }
    if (vertices_file.is_open()) {
      for (Vertex_index vert: vertexList) {
        Point_3 vPoint = mesh.point(vert);
        vertices_file << "{" << vPoint.x() << ", " << vPoint.y() << ", " << vPoint.z() << "}" << "\n";
      }
      vertices_file << "\n";
    }
    if (cmove_file.is_open()) {
      cmove_file << "{" << current_move.x() << ", " << current_move.y() << ", " << current_move.z() << "}";
      cmove_file << "\n";
    }

    std::pair<Point_3, std::vector<Vertex_index>> intersection_info = find_intersection(mesh, currentSourceFace, source_point, source_point+current_move, vertexList);
    Point_3 intersection_point = intersection_info.first;
    std::vector<Vertex_index> intersected_elements = intersection_info.second;
    
    if (intersection_point == Point_3(1000,1000,1000)) {
      intersection = false;

      std::cout << "no intersection, breaking loop" << std::endl;
      break;
    }

    Vector_3 vector_to_intersection = Vector_3(source_point, intersection_point);
    double lengthToSharedElement = vectorMagnitude(vector_to_intersection); //how far we've traveled
    current_move = reduceVector(current_move, lengthToSharedElement);         //decrease move size by length to intersected vertex/edge --
                                                                            //effectively the step where we "walk" to that intersection
    
    target = intersection_point+current_move; //storage of where the move vector currently points for rotation later

    currentTargetFace = getTargetFace(intersected_elements, source_point, vector_to_intersection, currentSourceFace, mesh); //face we're about to walk into;

    source_point = intersection_point;//update source to be the most recent intersection point -- finish walking there

    Face_location newMoveLocation = rotateIntoNewFace(mesh, currentSourceFace, currentTargetFace, source_point, target);
    Point_3 rotatedTarget = PMP::construct_point(newMoveLocation,mesh);
    std::cout << "rotated target to " << rotatedTarget << std::endl;

    //check that we've rotated in the right direction via overlap
    current_move = Vector_3(source_point, rotatedTarget);//source is now intersection

    overlap = current_move*currentTargetFaceNormal;
    if (overlap > baseAgreement+.01) {
    //PLACEHOLDER to see if sign convention is breaking
      std::cout << "reversing rotation angle!" << std::endl;
      rotatedTarget = rotateAboutAxis(forRotation, sharedEdge, -rotationAngle)[0];
      current_move = Vector_3(source_point, rotatedTarget);
    }
    overlap = current_move*currentTargetFaceNormal;
    if (overlap > baseAgreement +.01) {
      //unit test to see if we've got the shift right
      std::cout << "overlap too high to continue. " << std::endl;
      std::cout << "final target location " << rotatedTarget <<std::endl;
      break;
    }
    std::cout << "final target location " << rotatedTarget << std::endl;

    //check length of current_move before and after rotation
    std::cout << "overlap of current_move with target face normal: " << overlap  << std::endl;

    currentSourceFace = currentTargetFace;
  }
  //source_point+move is the location in the original face if there were no intersections, and it will 
  //be the location in the unfolded mesh if there were intersections (from an edge intersection to a spot
  //within a face)
  
  //declare final point & write it to file -- there is always one more path point than move vectors 
  Point_3 fTarget = source_point+current_move;
  if (intersections_file.is_open()){
    intersections_file << "{" << fTarget.x() << ", " << fTarget.y() << ", " << fTarget.z() << "}";
    intersections_file << "\n";
  }
  intersections_file.close();
  vertices_file.close();
  cmove_file.close();
  return fTarget; //might need some locating/more closely tying this to the mesh, but this should in principle be correct 
}
      
