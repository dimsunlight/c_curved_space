//Author: Toler H. Webb
//Description: code which takes as input a position on a mesh and a 
//displacement and returns as output a new position on the mesh. 
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
#include <math.h>
//types and names
typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::FT                                                           FT;
typedef K::Point_2                                                      Point_2;
typedef K::Ray_2                                                        Ray_2;
typedef K::Point_3                                                      Point_3;
typedef K::Vector_3                                                     Vector_3;
typedef K::Ray_3                                                        Ray_3;
typedef CGAL::Surface_mesh<Point_3>                                     Triangle_mesh;
typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor  vertex_descriptor;
typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor    face_descriptor;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;

//utility functions
auto normalize(Vector_3 v)
{
  auto const slen = v.x()*v.x() + v.y()*v.y()+v.z()*v.z();
  auto const d = CGAL::approximate_sqrt(slen);
  return v / d;
}

//function definitions
auto getVertexPositions(Triangle_mesh mesh, Triangle_mesh::face_index fIndex) {
   
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

auto getSharedEdge(std::vector<Point_3> vs1, std::vector<Point_3> vs2) {
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
  
  return sharedElements;
}

auto getVertexToRotate(std::vector<Point_3> vs2, std::vector<Point_3> sEdge) {
  //assuming second face includes the vertex we want to rotate
  std::vector<Point_3> oddVertex;
  
  for (Point_3 i: vs2) {
    if (std::find(sEdge.begin(), sEdge.end(), i) != sEdge.end()) {}
    else {
      oddVertex.push_back(i);
    }
  }

  return oddVertex;
}

auto rotateAboutSharedAxis(Point_3  target, std::vector<Point_3> axis, double angle) {
  //assumes that the shared axis is represented by a shared edge, rather than a vertex
  Vector_3 axisVector = normalize(Vector_3(axis[0],axis[1]));
  Point_3 shiftedTarget = target - Vector_3({0,0,0}, axis[0]);//can subtract vectors from points, but not points themselves
  std::cout << "target shifted for rotation to " << shiftedTarget << std::endl;
  /*
  double Identity[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  double crossMatrix[3][3] = {{},{},{}};
  double tensorProduct[3][3] = {{},{},{}};

  double formula = 1;
  */	
  return 0;
}

auto overEdge(Triangle_mesh mesh, Face_location f1, Face_location f2, Point_3 pos, Vector_3 move) {
  Face_location oldPosLocation = PMP::locate(pos,mesh);
  //compute normals are unit normals, so we shouldn't need additional normalization
  Vector_3 oldNormal = PMP::compute_face_normal(f1.first,mesh);
  Vector_3 newNormal = PMP::compute_face_normal(f2.first,mesh);
  //normals compute fine (verified via print) 
  std::cout << "Original face normal is " << oldNormal << std::endl;//to test if shift is orthogonal

  FT oldDotNew = CGAL::scalar_product(oldNormal,newNormal);
  double angle = acos(oldDotNew);
  //angle within expectations unit test 
  if (angle > 0.5) {
    std::cout << "Angle over 0.5!" << std::endl;
  }
  std::cout << "Angle: " << angle << std::endl;
  std::cout << "face 1: " << f1.first << std::endl;
  std::cout << "face 2: " << f2.first << std::endl;
  //task deconstruction for unfolding --essentially very simple version of geodesic finder
  std::vector<Point_3> vertices1 = getVertexPositions(mesh, f1.first);
  std::vector<Point_3> vertices2 = getVertexPositions(mesh, f2.first);
  std::cout << "face 1 vertices: " << std::endl; 
  for (Point_3 i: vertices1) std::cout << i << std::endl;
  std::cout << "face 2 vertices: " << std::endl;
  for (Point_3 i: vertices2) std::cout << i << std::endl;
  std::vector<Point_3> sharedEdge = getSharedEdge(vertices1,vertices2); 
  //by using Point_3 objects i'm opening myself up to numerical precision errors
  std::cout << "vertices of shared edge: " << std::endl;
  for (Point_3 i: sharedEdge) std::cout << i << std::endl;
  std::vector<Point_3> vertexToRotate = getVertexToRotate(vertices2, sharedEdge);
  std::cout << "target vertex coordinates to rotate" << std::endl;
  Point_3 vtr;
  for (Point_3 i: vertexToRotate){
	  std::cout << i << std::endl;
	  std::cout << typeid(i).name() << std::endl;
	  vtr = i;
  }
  //Point_3 newVertexLocation = rotateAboutSharedAxis(vtr, sharedEdge, angle);
  //tempFace = Face(f2.v1,f2.v2, newVertexLocation);
  //baryCoordinatesInNewFace = getBaryCoordinates(tempFace,pos+move); //coordinates tempFace are the same as those in f2
  //newXYZCoordinates = getXYZCoordinates(mesh,baryCoordinatesInNewFace,f2);
  //and done. 
  
  return Point_3(0.0,0.0,0.0);
}

auto shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
  //first: check if we will be in our out of triangle after naive shift
  //then: do the simple shift if we will, unfold shift if we won't 

  Point_3 trialNewPos = pos + move;
  Face_location oldPosLocation = PMP::locate(pos, mesh);
  Face_location newPosLocation = PMP::locate(trialNewPos,mesh);

  bool withinTriangle= oldPosLocation.first == newPosLocation.first;
  
  if (withinTriangle) {
    //could also return bary coordinates; PMP::locate(trialNewPos,mesh)
    std::cout << trialNewPos << std::endl;
    return trialNewPos;
  }

  if (!withinTriangle) {
    return overEdge(mesh, oldPosLocation, newPosLocation, pos, move);
  }

  std::cout << "neither in or out of triangle...?" << std::endl;
  return Point_3(0.0,0.0,0.0);
}


//main
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

 Vector_3 forceDisplacement = -0.1*Vector_3(-0.0145553, 0.0453073, 0.176461); //reversing direction because it looks more promising for tests when plotting everything in mathematica 
 Point_3 pointToMove = Point_3(3.51033, 1.9177, 0);
 


 Point_3 newPos = shift(tmesh, pointToMove, forceDisplacement);

 return 0;
}





