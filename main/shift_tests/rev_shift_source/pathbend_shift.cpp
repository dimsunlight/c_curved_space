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



face_descriptor getTargetFace(Segment_3 edge, Point_3 pos, Vector_3 toIntersection, face_descriptor currentSourceFace, Triangle_mesh mesh, bool throughVertex) {
 //goal is to find the face we're moving into. this is very easy if we share an edge, and harder if we only share a vertex. 
 //easy if we share an edge because we can just look for faces which share those vertices;
 //hard if we only share a vertex because the faces that share a vertex are not limited to the source and target.
 //maybe we can determine if we've gone through a vertex in advance. 
 if (throughVertex) {
   std::cout << "picking face -- through vertex..." << std::endl;

 }
 //big search function will be too slow... maybe?
 double moveEpsilon = 1.3; //using a tiny movement in the direction of the intersection vector to determine which face we're moving into

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

/*
Point_3 check2DIntersection(std::vector<Point_3> vertices, Segment_3 check segment) {
  
  Point_3 segment_source = check_segment.source();
  Vector_3 segment_vector = check_segment.to_vector();
  
  


  const auto result = CGAL::intersection(edge_segment,check_segment_2d);
  
  const Point_2* intersection_point = boost::get<Point_2>(&*result);
  if (result) {
    if(intersection_point) {
      std::cout << "INTERSECTION FOUND IN 2D" << std::endl;
      //now -- find distance to intersection, along the edge (which will be the 
      //first segment), and use the distance to find the equivalent point along the 
      //edge in the non-rotated space. 
      Vector_2 to_intersection = Vector_2(xys1s, *intersection_point);
      double lengthToIntersection = vectorMagnitude(to_intersection);
      Vector_3 vector_to_intersection = lengthToIntersection*normalizer(vector1); 
      std::cout << "Intersection Point: " << s1Source+vector_to_intersection << std::endl;
      return s1Source+vector_to_intersection;
    } else {
      std::cout << "line intersection" << std::endl;
      return Point_3(0,0,0);
    }
  } else {
    std::cout << "no intersection." << std::endl;
    return Point_3(0,0,0);
  }
}
*/

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
  bool throughVertexFlag;
  bool skip; //set this so we can ignore the remainder of a while loop once we've found an intersection

  //useful items for loop w/definition
  bool intersection = true; // true until we have checked all the edges/vertex and verified there's no intersection
  Segment_3 checkSegment = Segment_3(source_point, source_point+current_move); //check for edge intersections with this segment 
  std::cout << "Initial checksegment is " << checkSegment <<std::endl;

  while(intersection){
    vertexList = getVertexPositions(mesh,currentSourceFace);
    edgesList = createEdgeSegments(vertexList);
    skip = false; //this will continue to be redefined true as long as there is an intersection
    if (throughVertexFlag) throughVertexFlag = false;
    
    for (Segment_3 edge: edgesList) {
      
      std::cout << "searching for intersection with edge " << edge << std::endl;
      const auto result = CGAL::intersection(checkSegment,edge);    

      if (result) {
	
	const Point_3* edge_intersection = boost::get<Point_3>(&*result);
        	
        for (Point_3 vert: vertexList) {
	  if(CGAL::intersection(checkSegment,vert)) {
	    throughVertexFlag=true;
	    std::cout << "going through a vertex! hold on!" << std::endl;
	    continue;
	  }
	}
        std::cout << "Intersection, bending path" << std::endl;
	
	Vector_3 vector_to_intersection = Vector_3(source_point, *edge_intersection);
        double lengthToSharedElement = vectorMagnitude(vector_to_intersection); //how far we've traveled
        current_move = reduceMove(current_move, lengthToSharedElement); //decrease move size by length to intersected vertex/edge --
						       //effectively the step where we "walk" to that intersection

        source_point = *edge_intersection;//update source to be the most recent intersection point, as we have finished walking there
                                          //maybe need to pmp::locate -->pmp::construct point to make sure we stay on the mesh through numerical errors. Should be a tiny approximation at worst. 

        currentTargetFace = getTargetFace(edge, source_point, current_move, currentSourceFace, mesh, throughVertexFlag); //face we're about to walk into;
															 //works much better for smoother meshes
        currentTargetFaceNormal = PMP::compute_face_normal(currentTargetFace,mesh);
	std::cout << "bending the path" << std::endl;
    
        current_move = projectMoveDown(source_point, currentTargetFaceNormal, normalizer(current_move), vectorMagnitude(current_move)); //bend the path we are about to walk into the plane of the current face
        std::cout<< "Current move is " << current_move << std::endl;	
        checkSegment = Segment_3(source_point,source_point+current_move); //look ahead to see if we will intersect anything when we walk

        skip = true;	
        continue; // skip the rest of the for loop once we've found an intersection, if possible
      }
    } 
    if (skip) continue; //just avoiding last statement
    std::cout << "No intersection found. Shifting as if flat." << std::endl;
    intersection = false; //if we haven't found an intersection, there is no intersection. 
  }
  //source_point+move is the location in the original face if there were no intersections, and it will 
  //be the location in the unfolded mesh if there were intersections (from an edge intersection to a spot
  //within a face)  
  return source_point+current_move; //might need some locating/more closely tying this to the mesh, but this should in principle be completely correct 

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

 Vector_3 forceDisplacement = -2*Vector_3(-0.037457, 0.0330185, 0.177704); //reversing direction because it looks more promising for tests when plotting everything in mathematica 
 Point_3 pointToMove = Point_3(3.51033, 1.9177, 0);


 Point_3 newPos = shift(tmesh, pointToMove, forceDisplacement);
 std::cout << "The post-shift XYZ location is " << newPos << std::endl;
 return 0;
}


