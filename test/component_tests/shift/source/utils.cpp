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
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;

//utility functions
double vectorMagnitude(Vector_3 v) {
  auto const slen = v.x()*v.x() + v.y()*v.y()+v.z()*v.z();
  auto d = CGAL::approximate_sqrt(slen);
  return d; 
}

double vectorMagnitude(Vector_2 v) {
  auto const slen = v.x()*v.x() + v.y()*v.y();
  auto d = CGAL::approximate_sqrt(slen);
  return d;
}

Vector_3 normalizer(Vector_3 v)
{
  auto const d = vectorMagnitude(v);
  return v / d;
}

Vector_2 normalizer(Vector_2 v) 
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
std::vector<Point_3> getVertexPositions(Triangle_mesh mesh, Face_index fIndex) {

  Triangle_mesh::Halfedge_index hf = mesh.halfedge(fIndex);
  std::vector<Point_3> vertices;

  for(Triangle_mesh::Halfedge_index hi : halfedges_around_face(hf, mesh))
  {
    Triangle_mesh::Vertex_index vi = source(hi, mesh);
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
    Vertex_index vi = source(hi, mesh);
    vertices.push_back(vi);
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

double angleBetween(Face_index f1, Face_index f2, Triangle_mesh mesh) {
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

double angleBetween(Vector_3 v1, Vector_3 v2) { 
    v1 = normalizer(v1);
    v2 = normalizer(v2);
    double oneDotTwo = CGAL::scalar_product(v1,v2); 
    return acos(oneDotTwo); 
}

std::vector<double> baryCoords(Point_3 a, Point_3 b, Point_3 c, Point_3 p) {
  //use cramer's rule to determine the barycentric coordinates of p within the face described by a, b, c.
  //can also just use PMP::barycentric_coordinates, but I like having this on file.  
  Vector_3 v0 = Vector_3(a,b), v1 = Vector_3(a,c), v2 = Vector_3(a,p);
  //for cgal vector_3 objects, * indicates scalar product
  double d00 = v0*v0;
  double d01 = v0*v1;
  double d11 = v1*v1;
  double d20 = v2*v0;
  double d21 = v2*v1;
  double denom = d00*d01-d01*d01;
  double beta1 = (d11*d20-d01*d21)/denom;
  double beta2 = (d00*d21-d01*d20)/denom;
  double beta3 = 1.0-beta1-beta2;
  std::vector<double> baryCoords = {beta1,beta2,beta3};
  return baryCoords;
}

std::vector<Point_3> rotateAboutAxis(std::vector<Point_3> targets, std::vector<Point_3> axis, double rotAngle) {
  //rotation about an axis named axis; handles rotation of multiple targets at once
  
  Vector_3 axisVector = normalizer(Vector_3(axis[0],axis[1]));
  std::vector<Point_3> shiftedTargets; 
  shiftedTargets.reserve(3);

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
  rotatedTargets.reserve(3);
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

Point_3 rotateAboutAxis(Point_3 target, std::vector<Point_3> axis, double rotAngle) {
  //rotation about an axis named axis; handles rotation of multiple targets at once

  Vector_3 axisVector = normalizer(Vector_3(axis[0],axis[1]));
  Point_3 shiftedTarget = target-Vector_3({0,0,0},axis[0]);

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

  Point_3 rotatedTarget;
  double rotatedStorage[3];
  std::vector<double> targetStorage;

  //rotate and shift back all targets
  targetStorage = {shiftedTarget.x(),shiftedTarget.y(),shiftedTarget.z()};

  for (std::size_t i = 0; i < 3; i++) {
    double product = 0;
    for (std::size_t j = 0; j < 3; j++) {
      product += rMatrix[i][j]*targetStorage[j];
    }
    rotatedStorage[i] = product;
  }

  rotatedTarget = Point_3(rotatedStorage[0],rotatedStorage[1],rotatedStorage[2]);

  //undo shift at the end
  rotatedTarget = rotatedTarget + Vector_3({0,0,0}, axis[0]);

  return rotatedTarget;
}


auto createTemporaryMesh(std::vector<Point_3> vertexTrio) {
  Triangle_mesh tempMesh;
  Vertex_index u = tempMesh.add_vertex(vertexTrio[0]);
  Vertex_index v = tempMesh.add_vertex(vertexTrio[1]);
  Vertex_index w = tempMesh.add_vertex(vertexTrio[2]);
  Face_index  tFace = tempMesh.add_face(u,v,w);
  std::cout << "temp face normal = " << PMP::compute_face_normal(tFace,tempMesh) << std::endl;
  return tempMesh;
}

Point_3 getNewXYZLocation(Point_3 flatLocation, Triangle_mesh originalMesh, Triangle_mesh tempMesh, Face_index f2) {
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



Face_index getTargetFace(Point_3 pos, Vector_3 toIntersection, Face_index currentSourceFace, Triangle_mesh mesh) {
 //goal is to find the face we're moving into. this is very easy if we share an edge, and harder if we only share a vertex. 
 //easy if we share an edge because we can just look for faces which share those vertices;
 //hard if we only share a vertex because the faces that share a vertex are not limited to the source and target.
 
 //big search function will be too slow... maybe?
 double moveEpsilon = 1.05; //using a tiny movement in the direction of the intersection vector to determine which face we're moving into

 return PMP::locate(pos+moveEpsilon*toIntersection,mesh).first;//this is really, genuinely, just an approximation so i can debug the rest. 
                                                               //But it should identify things just fine most of the time. 
}

Vector_3 reduceVector(Vector_3 moveVector, double reduceLengthBy) {
  //reduce a vector by a length reduceLengthBy (for purposes of shift function, mostly)
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
  //we should probably just use CGAL's face normal finder routine instead most of the time; 
  //I know this version is meant to work in the mesh-structure-agnostic regime that I've been
  //operating in so far
  Point_3 projectedQuery;
  Vector_3 side1, side2, normal, toQuery;
  side1 = Vector_3(vertices[0],vertices[1]);
  side2 = Vector_3(vertices[0],vertices[2]);
  normal = CGAL::cross_product(side1, side2); //don't need to worry about direction due to dot
  toQuery = Vector_3(vertices[0], query);
  projectedQuery = vertices[0] + (toQuery-CGAL::scalar_product(toQuery,normal)*normal);
  return projectedQuery;
}

