//Author: Toler H. Webb
//Description: code which takes as input a position on a mesh and a
//displacement and returns as output a new position on the mesh.
//Uses barycentric coordinates to improve the generality of 
//intersection routines. 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
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
typedef CGAL::Surface_mesh<Point_3>                                     Triangle_mesh;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;
typedef Triangle_mesh::Face_index                                       Face_index;
typedef Triangle_mesh::Halfedge_index                                   Halfedge_index;
typedef Triangle_mesh::Vertex_index                                     Vertex_index; 
typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor    face_descriptor;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;

auto vector_magnitude(Vector_3 v) {
  auto const slen = v.x()*v.x() + v.y()*v.y()+v.z()*v.z();
  auto d = CGAL::approximate_sqrt(slen);
  return d;
}

double triangle_area(Point_3 v1, Point_3 v2, Point_3 v3) {
  Vector_3 side1 = Vector_3(v1, v2);
  Vector_3 side2 = Vector_3(v1, v3);
  return vector_magnitude(CGAL::cross_product(side1,side2))/2; 
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

std::vector<Point_3> getVertexPositions(Triangle_mesh mesh, face_descriptor fIndex) {

  Triangle_mesh::Halfedge_index hf = mesh.halfedge(fIndex);
  std::vector<Point_3> vertices; //I think we can say 3 because we know it's a triangle mesh

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

std::pair<Point_3,std::vector<Vertex_index>> find_intersection_full(Triangle_mesh mesh, Face_index sourceFace, Point_3 source, Point_3 target,  std::vector<Vertex_index> vertexIndices) {
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
  double tol = 0.0001; //checking for equivalence to zero

  //the first value is the only one not checked against tol in the main minimum-value-finding routine, so. 
  //we check it now. 
  if (toIntersect < tol) toIntersect = 2;//just needs to be some # greater than 1. 

  //min function with a max of 1. if we don't find something less than 1, no intersection.  
  for(double val: intersection_values) {
    std::cout << val << std::endl;
    if (val < toIntersect and val > tol) toIntersect = val;
  }

  std::cout << "toIntersect is " << toIntersect << std::endl;
    
  std::vector<Vertex_index> fillerVector = {}; //for null return later
  
  if (toIntersect < 1 and toIntersect > tol) {
    Point_3 min_intersection;
    min_intersection = Point_3(b11,b12,b13);
    min_intersection = min_intersection+toIntersect*displacement;
    std::array<double,3> newBaryCoords = {min_intersection.x(),min_intersection.y(),min_intersection.z()};
    
    //off vertices tell us what vertices are *not* part of the intersected edge
    int vertexIndicator = 0;
    std::vector<int> offVertices = {}; //indices of vertices with bary coordinate zero at the intersection
    //would use more economical "erase" for the below, but removing entries will mess up the indices if we have to erase twice
    std::cout << "bary coordinates of intersection: " << std::endl;
    for (double b: newBaryCoords) {
      std::cout << b << std::endl;
      if (b < tol) {
        offVertices.push_back(vertexIndicator);
      }
      vertexIndicator+=1;
    }
    std::cout << std::endl;

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
			       	min_intersection[2]*Vector_3(faceVertices[2].x(),faceVertices[2].y(),faceVertices[2].z()); //construct point from definition of barycentric to xyz conversion
    */
    std::cout << "intersection point " << std::endl;
    std::cout << xyz_intersection << std::endl; 

    return std::make_pair(xyz_intersection,intersected); // this version returns the barycentric intersection point
  }
  //if we don't have an intersection -- FIX CONDITIONALS LATER IN THE CODE
  else {
    std::cout << "No intersection. Returning filler point.";
    return std::make_pair(Point_3(1000,1000,1000),fillerVector);
  }
  return std::make_pair(Point_3(10000,1000,1000),fillerVector); //default return
}

Point_3 find_intersection(Vector_3 source_bary, Vector_3 query_bary, std::vector<Point_3> vertexList) {
  //returns barycentric coordinates of an edge intersection, if there is one
  Vector_3 displacement = query_bary - source_bary;
  double b11 = source_bary[0];
  double b12 = source_bary[1];
  double b13 = source_bary[2];
  std::vector<double> intersection_values = {-b11/displacement[0], -b12/displacement[1], -b13/displacement[2]}; 
  //all entries should by definition be positive, unless the source started outside of the triangle. 
  double toIntersect = intersection_values[0];
  for(double val: intersection_values) {
    if (val < toIntersect and val > 0) toIntersect = val;
  }
  std::cout << "toIntersect is " << toIntersect << std::endl; 
  if (toIntersect < 0) {
    std::cout << "distance along displacement is " << toIntersect << ". Is your source point correct?" << std::endl;
  }
  if (toIntersect < 1) { 
    std::cout << "found intersection. " << std::endl;
    Point_3 min_intersection;
    min_intersection = Point_3(b11,b12,b13);
    min_intersection = min_intersection+toIntersect*displacement;
    std::cout << "barycentric intersection point " << min_intersection << std::endl;    
    Vector_3 xyz_intersection = min_intersection[0]*Vector_3(vertexList[0].x(),vertexList[0].y(),vertexList[0].z()) + min_intersection[1]*Vector_3(vertexList[1].x(),vertexList[1].y(),vertexList[1].z())+ min_intersection[2]*Vector_3(vertexList[2].x(),vertexList[2].y(),vertexList[2].z());
    std::cout << "xyz intersection point (in f_i) " << xyz_intersection << std::endl;
    std::cout << " " << std::endl;
    return min_intersection; // this version returns the barycentric intersection point 
  }
  else {
    std::cout << "No intersection. Returning filler point."; 
    return Point_3(10,10,10); 
  }
  return Point_3(10,10,10);
}


Vector_3  move(Point_3 source, Point_3 target, Face_location source_mesh_coords, std::vector<Point_3> faceVertices, Triangle_mesh tmesh ) {
  Point_3 sourceDown = project_to_face(faceVertices, source);   
  Point_3 query = project_to_face(faceVertices, target);
  std::array<double, 3> source_bary = PMP::barycentric_coordinates(faceVertices[0],faceVertices[1],faceVertices[2], sourceDown); 
  std::array<double, 3> query_bary = PMP::barycentric_coordinates(faceVertices[0],faceVertices[1],faceVertices[2], query);
  Vector_3 source_bary_vector = Vector_3(source_bary[0],source_bary[1],source_bary[2]);
  Vector_3 query_bary_vector  = Vector_3(query_bary[0], query_bary[1] , query_bary[2]);
  std::cout << "source in barycentric " << source_bary_vector  <<   std::endl;
  std::cout << "target in barycentric " << query_bary_vector   <<   std::endl;
  Vector_3 bary_displacement_vector = query_bary_vector - source_bary_vector;

  std::cout << "displacement vector " <<   bary_displacement_vector << std::endl;
  std::cout << "displacement vector magnitude " << vector_magnitude(bary_displacement_vector) << std::endl;

  Point_3 intersection_coordinates = find_intersection(source_bary_vector, query_bary_vector, faceVertices); //using point_3 return type bc c++ whines about array returns
  std::cout << "bary after leaving f_i " << intersection_coordinates << std::endl;
  std::array<double,3> intersect_bary = {intersection_coordinates.z(), intersection_coordinates.x(), intersection_coordinates.y()};
  std::cout << "intersect_bary is " << intersect_bary[0] << " " << intersect_bary[1] << " " << intersect_bary[2] << std::endl;
  Face_location intersect_on_mesh = std::make_pair(source_mesh_coords.first, intersect_bary);
  Point_3 r3_intersection_coordinates = PMP::construct_point(intersect_on_mesh, tmesh);

  std::cout << "r3 intersection coordinates " << r3_intersection_coordinates << std::endl;
  return query_bary_vector;
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
 std::cout << typeid(tmesh.faces()).name() << std::endl;
 Vector_3 forceDisplacement = -.2*Vector_3(-0.037457, 0.0330185, 0.177704); //we start right next to an edge and this direction takes us over it
 Point_3 pointToMove = Point_3(3.51033, 1.9177, 0);
  
 Face_location onMesh = PMP::locate(pointToMove, tmesh); 
 Point_3 movedOntoMesh = PMP::construct_point(onMesh, tmesh);

 std::vector<Point_3> faceVertices = getVertexPositions(tmesh, onMesh.first);
 std::vector<Vertex_index> faceIndices = getVertexIndices(tmesh, onMesh.first);
 std::cout << "face vertices: " << faceVertices[0] << " " << faceVertices[1] << " " << faceVertices[2] << std::endl;
 std::cout << "face vertex indices: " << faceIndices[0] << " " << faceIndices[1] << " " << faceIndices[2] << std::endl;
 std::cout << " " << std::endl;
 std::cout << "source in xyz " << movedOntoMesh << std::endl;
 std::cout << "initial target in xyz " << movedOntoMesh + forceDisplacement << std::endl;
 std::cout << " " << std::endl;

 //move(movedOntoMesh,movedOntoMesh+forceDisplacement, onMesh, faceVertices, tmesh); // deprecated intersection test 
 
 std::pair<Point_3,std::vector<Vertex_index>> intersection_coordinates = find_intersection_full(tmesh, onMesh.first, movedOntoMesh, movedOntoMesh+forceDisplacement, faceIndices);
 

 //Point_3 xyzIntersectPoint = find_intersection_baryroutine(movedOntoMesh,movedOntoMesh+forceDisplacement, onMesh, faceVertices, tmesh);
 std::cout << "outputted xyz intersect point from full routine " << intersection_coordinates.first << std::endl;
 std::cout << "outputted vertices of intersection from full routine ";
 for (Vertex_index vi: intersection_coordinates.second) std::cout << vi << " "; 
 std::cout << std::endl;

 return 0;
}
