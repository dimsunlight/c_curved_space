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
typedef K::Segment_3                                                    Segment_3;
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
  if (sharedElements.size() == 1) std::cout << "Only 1 vertex shared." << std::endl;
  if (sharedElements.size() == 0) std::cout << "No shared vertices." << std::endl;
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

double angleBetween(Face_location f1, Face_location f2, mesh) {
  //generate normals -- compute_face_normal gives unit normals by default
  Vector_3 oldNormal = PMP::compute_face_normal(f1.first,mesh), newNormal = PMP::compute_face_normal(f2.first,mesh);
  //normals compute fine (verified via print) 

  FT oldDotNew = CGAL::scalar_product(oldNormal,newNormal);
  double angle = acos(oldDotNew);
  //angle within expectations unit test 
  if (angle > 1) {
    std::cout << "Angle between normals over 1!" << std::endl;
  }

  return angle;
}

auto verticesForUnfoldedFace(std::vector<Point_3> sharedEdge, std::vector<Point_3> vertices2 ) {
  //get original index locations -- fact that we always use second face is important here
  int ind1 = findIndex(sharedEdge[0],vertices2);
  int ind2 = findIndex(sharedEdge[1],vertices2); 
  //then, index of rotated vertex will be whichever of these is not included, so we create 
  //a phony index array to pick out the right one -- essentially just looks at the indices
  //1 and 2 and identifies which number between 0, 1, and 2 is not included
  std::vector<int> possibleIndices = {0,1,2};
  possibleIndices[ind1] = -1;
  possibleIndices[ind2] = -1;
  int rIndex;
  for (int ind: possibleIndices) {
    if (ind != -1) {
      rIndex = ind;
    }  
  } 

  //finally, tempTrio uses the found original indices to rearrange vertices into original order 
  std::vector<Point_3> tempTrio = {0,0,0};
  tempTrio[ind1] = sharedEdge[0];
  tempTrio[ind2] = sharedEdge[1];
  tempTrio[rIndex] = newVertexLocation;

  return tempTrio;
}

auto createTemporaryMesh(std::vector<Point_3> vertexTrio) {
  Triangle_mesh tempMesh;
  vertex_descriptor u = tempMesh.add_vertex(tempTrio[0]);
  vertex_descriptor v = tempMesh.add_vertex(tempTrio[1]);
  vertex_descriptor w = tempMesh.add_vertex(tempTrio[2]);
  face_descriptor tFace = tempMesh.add_face(u,v,w);
  std::cout << "temp face normal = " << PMP::compute_face_normal(tFace,tempMesh) << std::endl;
  return tempMesh;
}

Point_3 getNewXYZLocation(Point_3 flatLocation, Triangle_mesh tempMesh, Face_location f2) {
  Face_location tempLocation= PMP::locate(pos+move,tempMesh);
  Barycentric_coordinates b_weights = tempLocation.second;
  Face_location frankenLocation = std::make_pair(f2.first,b_weights);
  Point_3 newXYZLocation = PMP::construct_point(frankenLocation, mesh);

}

auto overEdge(Triangle_mesh mesh, Face_location f1, Face_location f2, Point_3 pos, Vector_3 move) {
  //functionality to move a particle location when it must travel over a mesh edge
  
  //first: grab original location via locate for use later
  Face_location oldPosLocation = PMP::locate(pos,mesh);

  //get angle between faces we need to move between
  double faceAngle = angleBetween(f1,f2,mesh);
  
  //grab the positions of the vertices for the two faces
  std::vector<Point_3> vertices1 = getVertexPositions(mesh, f1.first), vertices2 = getVertexPositions(mesh, f2.first);
  
  std::vector<Point_3> sharedEdge = getSharedEdge(vertices1,vertices2); 
  
  std::vector<Point_3> vertexToRotate = getVertexToRotate(vertices2, sharedEdge);
  //vertexToRotate returns a vector (in the case of a shared vertex but not a shared edge), but i haven't 
  //implemented behavior for that -- so for now we just grab the ``first'' element instead of splitting up
  //which functions we use in the following
  Point_3 newVertexLocation = rotateAboutSharedAxis(vertexToRotate[0], sharedEdge, -angle);
  
  std::vector<Point_3> tempTrio = verticesForUnfoldedFace(sharedEdge, vertices2);

  Triangle_mesh tempMesh = createTemporaryMesh(tempTrio);

  return getNewXYZLocation(pos+move,tempMesh, f2);
}

Point_3 shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
  //first: check if we will be in our out of triangle after naive shift
  //then: do the simple shift if we will, unfold shift if we won't 

  Face_location oldPosLocation = PMP::locate(pos, mesh);
  Point_3 meshPos = PMP::construct_point(oldPosLocation, mesh);
  Point_3 trialNewPos = meshPos + move;
  Face_location newPosLocation = PMP::locate(trialNewPos,mesh);

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


/*PLANNING VERSION OF NEW SHIFT
 * Point_3 shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
 *   oldPosLocation = PMP::locate(pos, mesh); //original position is fine
 *   //big departure from original -- we can't say anything about the new position until we've drawn the ray of length (len(move)) and learned about whether it does or doesn't intersect with an edge
 *   double TravelLength = vectorMagnitude(move);
 *   //now we need to see if the line segment from the original position to the position + the vector move (guaranteed to be in the tangent plane of the face) intersects any edges.
 *   //it's important that we use CGAL's version of the intersection routine here, rather than my own algebraic one. In three dimensions, homebrewed approaches are insufficiently
 *   //robust -- small numerical errors in heading will be detected as missed intersections. However, we *could* rewrite this whole setup as a two-dimensional one,
 *   //as all the pieces are really two-dimensional, and then try to do it from there with Segment_2 objects...      
 *   //oldPosLocation.first is the face the particle starts in.   
 *
 *   //we will eventually draw all vertex/edge segments of current face and store in lists; 
 *   std::vector<Point_3>    vertexList;
 *   std::vector<Segment_3> edgesList;
 *
 *   //create unfolded mesh structure to store faces as we unfold them 
 *   Triangle_mesh unfoldedMesh = createTemporaryMesh(vertexList);
 *   Face_descriptor currentFace = unfoldedMesh[0]; //should be only face of unfoldedMesh for now
 *   bool intersection = true; // true until we have checked all the edges/vertex and verified there's no intersection
 *   Segment_3 checkSegment = Segment_3(pos, pos+move); //now we define segment to be checked against... neat redefinitions later
 *   double lengthToSharedElement;
 *
 *   while(intersection){
 *     vertexList = getVertices(oldPosLocation.first);
 *     edgesList  = createEdgeSegments(vertexList);
 *     
 *     for (Point_3 vert: vertexList) {
 *       if (intersects(checkSegment,vert) {
 *         connectedFace = sharedFaces(vertex, oldPosLocation.first, face_list);
 *         std::vector<Point_3> ssvu = singleSharedVertexUnfolding(face1);
 *         currentFace = unfoldedMesh.add_face(ssvu[0],ssvu[1],ssvu[2]) //need to make sure these are in the original order!!!;
 *         move = reduceMove(move,lengthToSharedElement); //decrease move size by length to intersected vertex/edge. idk how we're getting that yet
 *         checkSegment = Segment_3(intersection_point, intersection_point+move); 
 *         continue; // skip the rest of the for & while loop once we've found an intersection, if possible
 *       }
 *     }
 *
 *     for (Segment_3 edge: edgesList) {
 *       if (intersects(checkSegment,edge) {
 *         connectedFace = sharedFace(edge, oldPosLocation.first, original_mesh_faces);
 *         std::vector<Point_3> seu = sharedEdgeUnfolding(currentFace,connectedFace,mesh); //mostly does what overEdge currently does
 *         currentFace = unfoldedMesh.add_face(seu[0],seu[1],seu[2]) //need to make sure these are in the original order!!!;
 *         move = reduceMove(move,lengthToSharedElement); //decrease move size by length to intersected vertex/edge. idk how we're getting that yet
 *         checkSegment = Segment_3(intersection_point, intersection_point+move); 
 *         continue; // skip the rest of the for & while loop once we've found an intersection, if possible
 *       }
 *     } 
 *  }
 * }
 */ 



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


 Point_3 newPos = shift(tmesh, pointToMove, forceDisplacement);
 std::cout << "The post-shift XYZ location is " << newPos << std::endl;
 return 0;
}



