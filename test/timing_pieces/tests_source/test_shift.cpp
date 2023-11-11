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
typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor  vertex_descriptor;
typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor    face_descriptor;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;

//utility functions
auto normalizer(Vector_3 v)
{
  auto const slen = v.x()*v.x() + v.y()*v.y()+v.z()*v.z();
  auto const d = CGAL::approximate_sqrt(slen);
  return v / d;
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

auto rotateAboutSharedAxis(Point_3 target, std::vector<Point_3> axis, double rotAngle) {
  //assumes that the shared axis is represented by a shared edge, rather than a vertex
  Vector_3 axisVector = normalizer(Vector_3(axis[0],axis[1]));
  Point_3 shiftedTarget = target - Vector_3({0,0,0}, axis[0]);//can subtract vectors from points, but not points themselves
  std::cout << "target shifted for rotation to " << shiftedTarget << std::endl;
  double u1 = axisVector.x(), u2 = axisVector.y(), u3=axisVector.z();
  
  double c = cos(rotAngle), s = sin(rotAngle);
  double C = 1-c;
  //unfortunately, I believe I have to multiply element wise in C++
  double Identity[3][3] = {{c,0.0,0.0},{0.0,c,0.0},{0.0,0.0,c}};
  double crossMatrix[3][3] = {{0,-s*u3,s*u2},{s*u3,0,-s*u1},{-s*u2,s*u1,0}};
  double tensorProduct[3][3] = {{C*u1*u1,C*u1*u2,C*u1*u3},{C*u2*u1,C*u2*u2,C*u2*u3},{C*u3*u1,C*u3*u2,C*u3*u3}};
  
  double rMatrix[3][3];
  //this double loop to define the rotation matrix is probably a place where it'd be easy to speed things up
  //with more c++ knowledge
  for (std::size_t i = 0; i < 3; i++) {
    for (std::size_t j = 0; j < 3; j++) {
	    rMatrix[i][j] = Identity[i][j] + crossMatrix[i][j] + tensorProduct[i][j];
    }
  }
  
  double rotatedStorage[3];
  double targetStorage[3] = {shiftedTarget.x(),shiftedTarget.y(),shiftedTarget.z()};
  
  for (std::size_t i = 0; i < 3; i++) {
    double product = 0;
    for (std::size_t j = 0; j < 3; j++) {
      product += rMatrix[i][j]*targetStorage[j];
    }
    rotatedStorage[i] = product; 
  }
  Point_3 rotatedTarget = Point_3(rotatedStorage[0],rotatedStorage[1],rotatedStorage[2]);
  
  rotatedTarget = rotatedTarget + Vector_3({0,0,0}, axis[0]);

  return rotatedTarget;
}

auto overEdge(Triangle_mesh mesh, Face_location f1, Face_location f2, Point_3 pos, Vector_3 move) {
  Face_location oldPosLocation = PMP::locate(pos,mesh);
  //compute normals are unit normals, so we shouldn't need additional normalization
  Vector_3 oldNormal = PMP::compute_face_normal(f1.first,mesh);
  Vector_3 newNormal = PMP::compute_face_normal(f2.first,mesh);
  //normals compute fine (verified via print) 

  FT oldDotNew = CGAL::scalar_product(oldNormal,newNormal);
  double angle = acos(oldDotNew);
  //angle within expectations unit test 
  if (angle > 1) {
    std::cout << "Angle over 1!" << std::endl;
  }

  //task deconstruction for unfolding --essentially very simple version of geodesic finder
  std::vector<Point_3> vertices1 = getVertexPositions(mesh, f1.first);
  std::vector<Point_3> vertices2 = getVertexPositions(mesh, f2.first);

  std::vector<Point_3> sharedEdge = getSharedEdge(vertices1,vertices2); 
  //by using Point_3 objects i'm opening myself up to numerical precision errors
  std::vector<Point_3> vertexToRotate = getVertexToRotate(vertices2, sharedEdge);
 
  Point_3 newVertexLocation = rotateAboutSharedAxis(vertexToRotate[0], sharedEdge, -angle);
  
  //get original index locations -- fact that we always use second face is important here
  int ind1 = findIndex(sharedEdge[0],vertices2);
  int ind2 = findIndex(sharedEdge[1],vertices2); 
  //then, index of rotated vertex will be whichever of these is not included, so we create 
  //a phony index array to pick out the right one -- essentially just looks at the indices
  //1 and 2 and identifies which number between 0, 1, and 2 is not included
  int possibleIndices[3] = {0,1,2};
  possibleIndices[ind1] = -1;
  possibleIndices[ind2] = -1;
  int rIndex;
  for (int ind: possibleIndices) {
    if (ind != -1) {
      rIndex = ind;
    }  
  } 

  //finally, tempTrio uses the found original indices to rearrange vertices into original order 
  Point_3 tempTrio[3];
  tempTrio[ind1] = sharedEdge[0];
  tempTrio[ind2] = sharedEdge[1];
  tempTrio[rIndex] = newVertexLocation;

  Triangle_mesh tempMesh;
  //through above tricky index manipulation, we should be able to guarantee that indices match
  //this creates a new mesh of a single face; barycentric weights in that face should match
  //the bary coordinates of the same face in the original mesh if the vertex order is the same. 
  vertex_descriptor u = tempMesh.add_vertex(tempTrio[0]);
  vertex_descriptor v = tempMesh.add_vertex(tempTrio[1]);
  vertex_descriptor w = tempMesh.add_vertex(tempTrio[2]);
  face_descriptor tFace = tempMesh.add_face(u,v,w);
  std::cout << "temp face normal = " << PMP::compute_face_normal(tFace,tempMesh) << std::endl;

  Face_location tempLocation= PMP::locate(pos+move,tempMesh);
  Barycentric_coordinates b_weights = tempLocation.second;
  
  Face_location frankenLocation = std::make_pair(f2.first,b_weights);
  
  Point_3 newXYZLocation = PMP::construct_point(frankenLocation, mesh);
  //baryCoordinatesInNewFace = getBaryCoordinates(tempFace,pos+move); //coordinates tempFace are the same as those in f2
  //newXYZCoordinates = getXYZCoordinates(mesh,baryCoordinatesInNewFace,f2);
  //and done. 
  return newXYZLocation;
}

Point_3 shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
  //first: check if we will be in our out of triangle after naive shift
  //then: do the simple shift if we will, unfold shift if we won't 
  clock_t preLocate  = clock();

  Face_location oldPosLocation = PMP::locate(pos, mesh);
  Point_3 meshPos = PMP::construct_point(oldPosLocation, mesh);
  Point_3 trialNewPos = meshPos + move;
  Face_location newPosLocation = PMP::locate(trialNewPos,mesh);
  std::cout << "Calls to locate took " << ((double) (clock() - preLocate)) /CLOCKS_PER_SEC << " seconds." << std::endl;

  bool withinTriangle= oldPosLocation.first == newPosLocation.first;
  
  if (withinTriangle) {
    //could also return bary coordinates; PMP::locate(trialNewPos,mesh)
    return trialNewPos;
  }

  if (!withinTriangle) {
    return overEdge(mesh, oldPosLocation, newPosLocation, meshPos, move);
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

 Vector_3 forceDisplacement = -0.12*Vector_3(-0.037457, 0.0330185, 0.177704); //reversing direction because it looks more promising for tests when plotting everything in mathematica 
 Point_3 pointToMove = Point_3(3.51033, 1.9177, 0);

 clock_t shiftStart = clock();
 Point_3 newPos = shift(tmesh, pointToMove, forceDisplacement);
 std::cout << "Shift took  " << ((double) (clock() - shiftStart)) / CLOCKS_PER_SEC << " seconds." << std::endl;

 std::cout << "The post-shift XYZ location is " << newPos << std::endl;
 return 0;
}

