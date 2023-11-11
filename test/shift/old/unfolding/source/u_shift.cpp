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

auto normalizer(Vector_3 v)
{
  auto const d = vectorMagnitude(v);
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



auto sharedEdge(Triangle_mesh mesh, Triangle_mesh unfolded_mesh, face_descriptor f1, face_descriptor f2, Point_3 pos, Vector_3 move,
                std::vector<Point_3> vertices1, std::vector<Point_3> vertices2, std::vector<Point_3> shared) {
  //functionality to move a particle location when it must travel over a mesh edge
  
  //get angle between faces we need to move between
  double faceAngle = angleBetween(f1,f2,mesh);

  std::vector<Point_3> vertexToRotate = getUnsharedElements(vertices2, shared);

  int rotatedIndex = findIndex(vertexToRotate[0],vertices2);
  int edgeIndex1 = findIndex(shared[0],vertices2);
  int edgeIndex2 = findIndex(shared[1], vertices2);
  
  //rotateAboutAxis is expanded to handle multiple targets for rotation, so 
  //it returns a std::vector -- that will be length 1 here and we can just pull the first element
  Point_3 newVertexLocation = rotateAboutAxis(vertexToRotate, shared, -faceAngle)[0];
  

  std::vector<Point_3> tempTrio = {Point_3(0,0,0),Point_3(0,0,0), Point_3(0,0,0)};

  tempTrio[rotatedIndex] = newVertexLocation;
  tempTrio[edgeIndex1] = shared[0];
  tempTrio[edgeIndex2] = shared[1]; 
  
  vertex_descriptor v0 = unfolded_mesh.add_vertex(tempTrio[0]);
  vertex_descriptor v1 = unfolded_mesh.add_vertex(tempTrio[1]);
  vertex_descriptor v2 = unfolded_mesh.add_vertex(tempTrio[2]);
 
  unfolded_mesh.add_face(v0,v1,v2);
  
  return unfolded_mesh;
}

auto sharedVertex(Triangle_mesh mesh, Triangle_mesh unfolded_mesh, face_descriptor f1, face_descriptor f2, Point_3 pos, Vector_3 move,
                std::vector<Point_3> vertices1, std::vector<Point_3> vertices2, std::vector<Point_3> shared) {
  //functionality to unfold if we must travel through a vertex into a face that only shares that vertex 
  
  //get angle between faces we need to move between
  double faceAngle = angleBetween(f1,f2,mesh);
  
  std::vector<Point_3> verticesToRotate = getUnsharedElements(vertices2,shared); //two vertices to rotate from new face...
  
  std::vector<Point_3> unsharedEdge = getUnsharedElements(vertices1,shared); //are rotated around the axis of the two unshared vertices from the original. 
  
  int rotatedIndex1 = findIndex(verticesToRotate[0],vertices2);
  int rotatedIndex2 = findIndex(verticesToRotate[1],vertices2);
  int staticVertexIndex = findIndex(shared[0],vertices2);

  std::vector<Point_3> newVertexLocations = rotateAboutAxis(verticesToRotate, unsharedEdge, -faceAngle); //should this be -faceangle? 
  
  std::vector<Point_3> tempTrio = {Point_3(0,0,0),Point_3(0,0,0), Point_3(0,0,0)};

  tempTrio[rotatedIndex1] = newVertexLocations[0];
  tempTrio[rotatedIndex2] = newVertexLocations[1];
  tempTrio[staticVertexIndex] = shared[0]; 

  vertex_descriptor v0 = unfolded_mesh.add_vertex(tempTrio[0]);
  vertex_descriptor v1 = unfolded_mesh.add_vertex(tempTrio[1]);
  vertex_descriptor v2 = unfolded_mesh.add_vertex(tempTrio[2]);

  unfolded_mesh.add_face(v0,v1,v2);

  return unfolded_mesh;
}  

auto unfold(Triangle_mesh mesh, Triangle_mesh unfolded_mesh, face_descriptor f1, face_descriptor f2, Point_3 pos, Vector_3 move) {
  //folding function that calls either sharedEdge or sharedVertex depending on the output of getSharedElements.
  //Returns the new Point_3 location after moving into the appropriate new face.  
  
  //grab the positions of the vertices for the two faces
  std::vector<Point_3> vertices1 = getVertexPositions(mesh, f1);
  std::vector<Point_3> vertices2 = getVertexPositions(mesh, f2);
  
  std::vector<Point_3> sharedElements = getSharedElements(vertices1,vertices2); 
  
  if (sharedElements.size() == 2) {
    return sharedEdge(mesh, unfolded_mesh, f1, f2, pos, move, vertices1, vertices2, sharedElements); 
  }
  
  if (sharedElements.size() == 1) {
    return sharedVertex(mesh, unfolded_mesh, f1, f2, pos, move, vertices1, vertices2, sharedElements); 
  }
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
 double moveEpsilon = 1.1; //using a tiny movement in the direction of the intersection vector to determine which face we're moving into

 return PMP::locate(pos+moveEpsilon*toIntersection,mesh).first;//this is really, genuinely, just an approximation so i can debug the rest. 
}

Vector_3 reduceMove(Vector_3 moveVector, double reduceLengthBy) {
  double originalLength = vectorMagnitude(moveVector);
  Vector_3 normalizedMove = normalizer(moveVector); 
  Vector_3 reducedMove = (originalLength-reduceLengthBy)*normalizedMove;

  return reducedMove;
}


Point_3 shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
  //big departure from original -- we can't say anything about the new position until we've drawn the ray of length
  //(len(move)) and learned about whether it does or doesn't intersect with an edge
  double travelLength = vectorMagnitude(move);
   
  //we will eventually draw all vertex/edge segments of current face and store in lists; 
  Face_location sourceLocation = PMP::locate(pos, mesh);
  Point_3 source_point = PMP::construct_point(sourceLocation, mesh); 
  face_descriptor currentSourceFace = sourceLocation.first;

  //initializations
  face_descriptor currentTargetFace;
  double lengthToSharedElement;
  std::vector<face_descriptor> faceIndexList; //store face indices here so we know where to look later
  bool skip; //set this so we can ignore the remainder of a while loop once we've found an intersection

  std::vector<Point_3> vertexList;
  std::vector<Segment_3> edgesList; 

  //useful items for loop w/definition
  bool intersection = true; // true until we have checked all the edges/vertex and verified there's no intersection
  Segment_3 checkSegment = Segment_3(pos, pos+move); //check for edge intersections with this segment 
  vertexList = getVertexPositions(mesh,currentSourceFace);

  Triangle_mesh unfoldedMesh = createTemporaryMesh(vertexList);
 
  bool throughVertexFlag;

  std::size_t numberUnfoldings = 0;

  
  while(intersection){
    vertexList = getVertexPositions(mesh,currentSourceFace);
    edgesList = createEdgeSegments(vertexList);
    skip = false; //this will continue to be redefined true as long as there is an intersection
    if (throughVertexFlag) throughVertexFlag = false;
    faceIndexList.push_back(currentSourceFace); //last entry of faceIndexList is always true index of most 
                                                //recently unfolded face 
    
    for (Segment_3 edge: edgesList) {

      const auto result = CGAL::intersection(checkSegment,edge);    //can handle motion along an edge with a conditional later !!!!!!!!
							      //eg: ifLineIntersection(checkSegment,edge) can handle both moving along the
							      //edge but within its length and moving along the edge and past an endpoint 
      if (result) {
        if(const Segment_3* intersection_segment = boost::get<Segment_3>(&*result)) {
          std::cout << "SEGMENT LINE INTERSECTION, EXITING"  << std::endl;
	  break;
        }
	
	const Point_3* edge_intersection = boost::get<Point_3>(&*result);
       
        	
        for (Point_3 vert: vertexList) {
	  if(CGAL::intersection(checkSegment,vert)) {
	    throughVertexFlag=true;
	    std::cout << "going through a vertex! hold on!" << std::endl;
	    continue;
	  }
	}
        
	//SO CORE ISSUE
	//we are trying to walk both in the unfolded mesh and the regular mesh
	//we need to track our position in both 
	//at the moment, the procedure is pretty complex and difficult to untangle
	//tomorrow: try writing path folding (rather than unfolding) version.
	//saturday, start untangling
	//sunday, hopefully finish. 

	Vector_3 vector_to_intersection = Vector_3(source_point, *edge_intersection);
	
	double lengthToSharedElement = vectorMagnitude(vector_to_intersection);  

        currentTargetFace = getTargetFace(edge, source_point, move, currentSourceFace, mesh, throughVertexFlag); 
	//technically, this can only give a face index, and we have to construct the face location using the updatedPos
	std::cout << "performing unfolding" << std::endl;
	//unfold adds the next unfolding to the temporary mesh. We can then get the position of the intersection point in that. 
	unfoldedMesh = unfold(mesh, unfoldedMesh, currentSourceFace, currentTargetFace, pos, move);
        numberUnfoldings += 1;
    
        source_point = *edge_intersection;//update source to be the most recent intersection point
	move = reduceMove(move,lengthToSharedElement); //decrease move size by length to intersected vertex/edge
        checkSegment = Segment_3(source_point,source_point+move); //change checkSegment to be the segment from the intersection point forward
	currentSourceFace = currentTargetFace; 

        skip = true;	
        continue; // skip the rest of the for loop once we've found an intersection, if possible
      }
    } 
    if (skip) continue; //just avoiding last statement
    intersection = false; //if we haven't found an intersection, there is no intersection. 
  }
  //source_point+move is the location in the original face if there were no intersections, and it will 
  //be the location in the unfolded mesh if there were intersections (from an edge intersection to a spot
  //within a face)  
  return getNewXYZLocation(source_point+move, mesh, unfoldedMesh, faceIndexList.back()); 

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

 Vector_3 forceDisplacement = -0.12*Vector_3(-0.037457, 0.0330185, 0.177704); //reversing direction because it looks more promising for tests when plotting everything in mathematica 
 Point_3 pointToMove = Point_3(3.51033, 1.9177, 0);


 Point_3 newPos = shift(tmesh, pointToMove, forceDisplacement);
 std::cout << "The post-shift XYZ location is " << newPos << std::endl;
 return 0;
}


/*
auto verticesForUnfoldedFace(std::vector<Point_3> shared, std::vector<Point_3> vertices2, std::vector<Point_3> newLocations ) {
  //defunct function for getting vertices in original order 
  std::vector<int> indices; 
  for (Point_3 sharedVertex: shared) {
    indices.push_back(findIndex(sharedVertex,vertices2); 
    //first element of shared will have index indices[0], second will have index indices[1], etc. 
  }

  //to handle case where we've rotated w/only one shared element, we can't use process of elimination
  std::vector<int> possibleIndices = {0,1,2};
  for (int index: indices) possibleIndices[index] = -1; 

  int rIndex;
  for (int ind: possibleIndices) {
    if (ind != -1) {
      rIndex = ind;
    }  
  } 

  //finally, tempTrio uses the found original indices to rearrange vertices into original order 
  std::vector<Point_3> tempTrio = {Point_3(0,0,0),Point_3(0,0,0),Point_3(0,0,0)};//null points for init. so we can have size 3
  tempTrio[ind1] = sharedEdge[0];
  tempTrio[ind2] = sharedEdge[1];
  tempTrio[rIndex] = newLocations[0];

  return tempTrio;
}
*/


/*
Point_3 shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move) {
  Face_location oldPosLocation = PMP::locate(pos, mesh); //original position is fine
  //big departure from original -- we can't say anything about the new position until we've drawn the ray of length
  //(len(move)) and learned about whether it does or doesn't intersect with an edge
  double TravelLength = vectorMagnitude(move);
   
  //oldPosLocation.first is the face the particle starts in.   

  //we will eventually draw all vertex/edge segments of current face and store in lists; 
  std::vector<Point_3>   vertexList;
  std::vector<Segment_3> edgesList;

  vertexList = getVertexPositions(mesh)
  //initializations
  face_descriptor connectedFace;
  double lengthToSharedElement;
  bool skip; //set this so we can ignore the remainder of a while loop once we've found an intersection

  //useful items for loop w/definition
  bool intersection = true; // true until we have checked all the edges/vertex and verified there's no intersection
  Segment_3 checkSegment = Segment_3(pos, pos+move); //now we define segment to be checked against... redefine later
  Triangle_mesh unfoldedMesh = createTemporaryMesh(vertexList);
  Face_descriptor currentFace = unfoldedMesh[0]; //should be only face of unfoldedMesh for now

  while(intersection){
    std::vector<Segment_3> edgesList = createEdgeSegments(vertexList);
    skip = false;

    for (Segment_3 edge: edgesList) {
      auto result = CGAL::intersection(checkSegment,edge);

      if (intersects(checkSegment,edge)) {
        connectedFace = sharedFace(edge, oldPosLocation.first, original_mesh_faces);
        std::vector<Point_3> seu = sharedEdgeUnfolding(currentFace,connectedFace,mesh); //mostly does what overEdge currently does
        currentFace = unfoldedMesh.add_face(seu[0],seu[1],seu[2]); //need to make sure these are in the original order!!!;
        move = reduceMove(move,lengthToSharedElement); //decrease move size by length to intersected vertex/edge. idk how we're getting that yet
        checkSegment = Segment_3(intersection_point, intersection_point+move);
        skip = true;	
        continue; // skip the rest of the for & while loop once we've found an intersection, if possible
      }
    } 
    if (skip) continue; //just avoiding last statement
    intersection = false; //if we haven't found an intersection, there is no intersection. 
  }
}
*/ 
