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
#include <CGAL/boost/graph/iterator.h>
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
typedef CGAL::Face_around_target_circulator<Triangle_mesh>              Face_circulator;
typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor           vertex_descriptor;
typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor             face_descriptor;

std::pair<Point_3,std::vector<Vertex_index>> find_intersection(Triangle_mesh mesh, Face_index sourceFace, Point_3 source, Point_3 target,  std::vector<Vertex_index> vertexIndices) {
  //use the barycentric coordinates of a point in order to determine the R3 coordinates of its intersection with an edge of the current face, if such an intersection exists.
  //takes only the r3 positions of the source and prospective target along with the vertices of the current face; intersection test exists independent of the mesh.
  
  std::vector<Point_3> faceVertices;
  for (Vertex_index vi: vertexIndices) faceVertices.push_back(mesh.point(vi));
  std::cout << "original vertex indices for intersection routine: " << std::endl;
  for (Vertex_index vi: vertexIndices) std::cout << vi << ",";
  std::cout << std::endl;  
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
  double tol = pow(10,-6); //checking for equivalence to zero within reasonable error

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
    for (double val: newBaryCoords) {
      std::cout << val << ",";	
    }
    std::cout << std::endl; 
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

    for (int i=0; i < vertexIndices.size(); i++) {
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

Face_index selectFaceFromVertex(const Vertex_index &intersectedVertex, const Vector_3 &toIntersection, const Face_index &source_face,
	       const Triangle_mesh &mesh) {
  //key piece here is the "face_around_target" circulator from CGAL

  std::vector<Face_index> candidateFaces; 
  candidateFaces.reserve(6); 
  
  //get all possible faces that we could walk into. 
  Face_circulator fbegin(mesh.halfedge(intersectedVertex),mesh), done(fbegin); //fbegin is the name of the Face_circulator type object, things 
  										 //in the parentheses are declaration arguments for Face_circulators
  Face_index cand; 
  do {
      cand = *fbegin++; //the next face to add is the face pointed to by the next
      			//iteration of the iterator 
      if (cand == source_face) continue; //don't need to consider the face we're coming from
      else {
        candidateFaces.push_back(cand);
      }
    } while(fbegin != done);

  std::vector<Point_3> faceVertices;

  std::cout << "selecting face given intersected vertex! " << intersectedVertex << std::endl;

  Point_3 source_r3 = mesh.point(intersectedVertex);
  // small bandaid -- because pmp::locate_vertex only works with boost graph vertex_descriptors, 
  // and we're using Vertex_index objects for type consistency, I first get the exact r3 location 
  // of the vertex above and then find its mesh location using PMP::locate. If we had written
  // the whole code with vertex_descriptors I wouldn't have to do this, but that would mess up 
  // ranging in other places.   
  Face_location intersectedBary = PMP::locate(source_r3,mesh);  

  // now -- check every candidate face to see if the path falls within it after bending. 
  // hypothesis: the path will only fall within the candidate face if that face is the *correct* face, most of the time -- 
  // pathological meshes where two faces are at an extreme angle relative to each other & very thin can break this. 
  // But if that happens, we're probably getting very nearly best possible behavior anyway. See through_vertex.nb notebook
  // for simple tests of this hypothesis.
  // We can check by extending the path out and seeing if its endpoint either intersects the candidate face or falls within
  // that face. If either are true, this is the correct face and we can break the loop.
  Face_index correctFace = source_face;
  bool outOfTriangle;
  std::array<double, 3> source_bary = intersectedBary.second;
 
  // we could while loop instead of using break;, but for loop guaranteed terminates
  for (Face_index candidate_face: candidateFaces) {
    std::vector<Vertex_index> faceVertices = getVertexIndices(mesh, candidate_face); 
    Face_location candidateLocation = rotateIntoNewFace(mesh, source_face, candidate_face, source_r3, source_r3 +  toIntersection);
    
    std::array<double, 3> cand_bary = candidateLocation.second;

    outOfTriangle = false;
    for (int i = 0; i < 3; i++) { 
      if (cand_bary[i] < 0) outOfTriangle = true;
    } 
    
    if (!outOfTriangle) {
      correctFace = candidate_face;
      break;
    }

    if (outOfTriangle) {      
      Point_3 cand_r3 = PMP::construct_point(candidateLocation, mesh); 
      std::pair<Point_3,std::vector<Vertex_index>> ifIntersect = find_intersection(mesh, source_face, source_r3, cand_r3, faceVertices);
      Point_3 intersection_point = ifIntersect.first; 
      if (intersection_point != Point_3(1000,1000,1000)) { 
	correctFace = candidate_face;
	break; 
      }
    } 
  }
  
  // we should always find the right face here, given that we're staying on the mesh. 
  // If we don't find one, something's wrong, and I want to know
  if (correctFace == source_face) {
    std::cout << "Found correct face to be source face, DEBUG" << std::endl;
  }	
  std::cout << "Using throughvertex routine, found target face: " << correctFace << std::endl;
  return correctFace;
}

/**
 * Find the face shift is about to travel into as it wraps the path around the mesh. 
 * params: 
 * 	intersected: the list of intersected vertices. If two, represents an edge; if one, represents a single vertex
 *	toIntersection: the Vector_3 object pointing from pos to the intersection point. Goes in the same direction as move vector.
 *	source_face: the face we're coming from.
 *	mesh: the mesh we're in, for reference. Should pass by address.  
 * returns the target face if one exists.  
 * 	
 */

Face_index getTargetFace(std::vector<Vertex_index> intersected, const Vector_3 &toIntersection, const Face_index &source_face, const Triangle_mesh &mesh) {
 //Find the face we're moving into by evaluating the other face attached to the edge made of vertex indices "intersected."

  std::cout << "printing intersected size: " << intersected.size() << std::endl;
  bool throughVertex = (intersected.size() == 1);
  std::cout << "through vertex? " << throughVertex << std::endl;
  if (throughVertex) {
    Vertex_index intersectedVertex = intersected[0];
    // usually, below picks from six possible faces. 
    return selectFaceFromVertex(intersectedVertex, toIntersection, source_face, mesh);   
  }
  Halfedge_index intersected_edge = mesh.halfedge(intersected[0],intersected[1]);
  Face_index provisional_target_face = mesh.face(intersected_edge);
  if (provisional_target_face == source_face) provisional_target_face = mesh.face(mesh.opposite(intersected_edge));
  return provisional_target_face;
}


Point_3 shift(const Triangle_mesh &mesh, const Point_3 &pos, const Vector_3 &move) {
  double travelLength = vectorMagnitude(move);

  //short pseudo-projection routine in case we're slightly off the face (so code is general to arbitrary starting point) 
  Face_location sourceLocation = PMP::locate(pos, mesh);
  Point_3 source_point = PMP::construct_point(sourceLocation, mesh); //xyz point representing current source

  Face_index currentSourceFace = sourceLocation.first;  
  Vector_3 currentSourceNormal = PMP::compute_face_normal(currentSourceFace,mesh);
 
  Point_3 target = source_point+move;

  std::vector<Vertex_index> vertexList;
  std::vector<Point_3> vertexPos = getVertexPositions(mesh, sourceLocation.first);

  //if the movement vector isn't quite in the tangent plane, put it there. 
  if (abs(currentSourceNormal*move) > 0) {
    target = project_to_face(vertexPos, pos+move);
  }
  Vector_3 current_move = Vector_3(source_point, target); 
  
  //if there is no intersection, avoid initalizing any of the intersection code
  std::array<double,3> targetBary = PMP::barycentric_coordinates(vertexPos[0],vertexPos[1],vertexPos[2],target);
  bool intersection = false;
  
  for (double bc: targetBary) {  
    if (bc < 0) {
      intersection = true; 
    }
  }

  if (intersection == false) {
     return source_point + current_move; 
  }

  //initializations if there *is* an intersection
  std::vector<Point_3> targetVertices;
  std::vector<Point_3> forRotation;
  Point_3 rotatedTarget;
  Face_index currentTargetFace;
  //Vector_3 currentTargetFaceNormal;
  double lengthToSharedElement;
  double rotationAngle;
  double overlap;

  while(intersection){
    //vertexList = getVertexPositions(mesh,currentSourceFace);
    vertexList = getVertexIndices(mesh,currentSourceFace);


    std::pair<Point_3, std::vector<Vertex_index>> intersection_info = find_intersection(mesh, currentSourceFace, source_point, source_point+current_move, vertexList);
    Point_3 intersection_point = intersection_info.first;
    std::vector<Vertex_index> intersected_elements = intersection_info.second;
    
    if (intersection_point == Point_3(1000,1000,1000)) {
      intersection = false;
      break;
    }

    Vector_3 vector_to_intersection = Vector_3(source_point, intersection_point);
    double lengthToSharedElement = vectorMagnitude(vector_to_intersection); //how far we've traveled
    current_move = reduceVector(current_move, lengthToSharedElement);         //decrease move size by length to intersected vertex/edge --
                                                                            //effectively the step where we "walk" to that intersection
    
    target = intersection_point+current_move; //storage of where the move vector currently points for rotation later
    
    currentTargetFace = getTargetFace(intersected_elements, vector_to_intersection, currentSourceFace, mesh); //face we're about to walk into;
    
    source_point = intersection_point;//update source to be the most recent intersection point -- finish walking there

    Face_location newMoveLocation = rotateIntoNewFace(mesh, currentSourceFace, currentTargetFace, source_point, target);
    Point_3 rotatedTarget = PMP::construct_point(newMoveLocation,mesh);

    //check that we've rotated in the right direction via overlap
    current_move = Vector_3(source_point, rotatedTarget);//source is now intersection
    currentSourceFace = currentTargetFace;
  }
  //source_point+move is the location in the original face if there were no intersections, and it will 
  //be the location in the unfolded mesh if there were intersections (from an edge intersection to a spot
  //within a face)
  
  Point_3 fTarget = source_point+current_move;
  return fTarget;  
}


//unit test code that can be thrown back in 
    //overlap = current_move*currentTargetFaceNormal;
    //std::cout << "overlap is " << overlap << std::endl; //if your shifts get weird, uncomment this; if the 
    //numbers are significantly deviating from zero, fix the below test and try that to see if angle 
    //conventions are getting broken. 

    //below unit test can cause spurious failures if the base overlap is the worst (which it often is
    //when feeding in testing moves). sharedEdge doesn't exist at present -- can rewrite to use 
    //face normals to create a rotation axis. 
    /*
    if (overlap > baseAgreement+.001) {
    //PLACEHOLDER to see if sign convention is breaking
      std::cout << "overlap is " << overlap << std::endl;
      std::cout << "reversing rotation angle!" << std::endl;
      //in below, sharededge doesn't exist -- this is actually just  a debug message
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
    */


