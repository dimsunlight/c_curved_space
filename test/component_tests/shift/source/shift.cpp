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
typedef CGAL::Vertex_around_target_circulator<Triangle_mesh>            Vertex_circulator;
typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor           vertex_descriptor;
typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor             face_descriptor;


std::pair<Point_3,std::vector<Vertex_index>> find_intersection(Triangle_mesh mesh, Face_index sourceFace, Point_3 source, Point_3 target,  std::vector<Vertex_index> vertexIndices) {
  //use the barycentric coordinates of a point in order to determine the R3 coordinates of its intersection with an edge of the current face, if such an intersection exists.
  //takes only the r3 positions of the source and prospective target along with the vertices of the current face; intersection test exists independent of the mesh.
  
  double tol = pow(10,-16); //checking for equivalence to zero within reasonable error
  //what do we do if we're heading onto another face but we're currently on an edge? We only want to log 
  //the intersection if we're specifically on the edge that is shared with the new face, 
  //but we do need to log that intersection automatically... 

  std::vector<Point_3> faceVertices;
  for (Vertex_index vi: vertexIndices) faceVertices.push_back(mesh.point(vi));
  //std::cout << "original vertex indices for intersection routine: ";
  //for (Vertex_index vi: vertexIndices) std::cout << vi << ",";
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
  
  //the first value is the only one not checked against tol in the main minimum-value-finding routine, so. 
  //we check it now. 
  if (toIntersect < -tol) toIntersect = 2;//just needs to be some # greater than 1. 

  //min function with a max of 1. if we don't find something less than 1, no intersection. 
  printf("printing vals from toIntersect:\n");
  for(double val: intersection_values) {
    std::cout << val << std::endl; 
    if (val < toIntersect and val > -tol) toIntersect = val;
  }  
  std::vector<Vertex_index> fillerVector = {}; //for null return later
 

  if (toIntersect < 1 and toIntersect > tol) {
    Point_3 min_intersection;
    min_intersection = Point_3(b11,b12,b13);
    min_intersection = min_intersection+toIntersect*displacement;
    
    // new bary coords indicates a trio of solutions to linear equations where one of the barycentric coords 
    // should be zero along the set displacement "displacement". Using point_3 object for convenience
    std::array<double,3> newBaryCoords = {min_intersection.x(),min_intersection.y(),min_intersection.z()};
    //print bary coordinates for check
    // off vertices tell us what vertices are *not* part of the intersected edge -- they will have value 0 
    int vertexIndicator = 0;
    std::vector<int> offVertices = {}; // indices of vertices with bary coordinate zero at the intersection
    // would use more economical "erase" for the below, but removing entries will mess up the indices 
    // if we have to erase twice, and we don't know a priori which index is being removed/in what order. 
    
    // below logic: We figure out which vertices have zero barycentric weight, then exclude them 
    // from the list of intersected vertices, which will collectively define an edge (two nonzero) 
    // or vertex (one nonzero) object. Way of working around index and position info not being
    // in the same object, and not stored the same way in this fxn 

    // if the barycentric coordinate that we found using minimal intersection is 0, that means 
    // it is not part of the intersected object (edge or vertex). Below finds which bary coords 
    // are zero so we can exclude them later. vertexIndicator is a standard counter so I can use
    // comprehension for the loop -- could have just as easily done a (for(int ii = 0)... ) style 
    // loop instead
    
    std::cout << "Barycentric intersection coordinates: "; 
    for (double b: newBaryCoords) {
      std::cout << b << ", "; 
      if (b < tol) {
        offVertices.push_back(vertexIndicator);
      }
      vertexIndicator+=1;
    }
    printf("\n"); 

    // now -- for the vertices that did *not* have barycentric weights = 0, we need to 
    // fill the intersected vector 
    std::vector<Vertex_index> intersected = {};
    bool flag = false;

    // loop through all the vertex indices 
    for (int i=0; i < vertexIndices.size(); i++) {
      // search the list of vertices which are not intersected (via 0,1,2 array indices) 
      // and see if this is one of them
      for (int off: offVertices) {
        if (i==off) flag = true;
      }

      // if it was one of the vertices that wasn't intersected (which is what we 
      // actually track), reset the flag and move to the next iteration
      if (flag == true) {
        flag = false;
	continue;
      }
      // if this vertex was intersected (nonzero bary weight), add it to the list
      // of intersected objects 
      intersected.push_back(vertexIndices[i]);
    }

    Face_location newPosition = std::make_pair(sourceFace,newBaryCoords);
    Point_3 xyz_intersection = PMP::construct_point(newPosition,mesh);
    return std::make_pair(xyz_intersection,intersected); // this version returns the xyz intersection point
  }
  //if we don't have an intersection, return filler point; for later, CGAL i think also has special Point_3 null types
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

void clampBary(const Triangle_mesh &mesh, const Face_index &face, Barycentric_coordinates &baryCoords) {
  // modify a set of bary coordinates in place by clamping them into a face. 	
  std::vector<Point_3> vList = getVertexPositions(mesh, face); 
  Face_location floc = std::make_pair(face, baryCoords);
  Vector_3 f_normal = PMP::compute_face_normal(face,mesh);
  Vector_3 toQuery = Vector_3(vList[0], PMP::construct_point(floc,mesh));
  Point_3 projectedQuery = vList[0] + (toQuery - CGAL::scalar_product(toQuery,f_normal)*f_normal);
  Barycentric_coordinates projected_barys = PMP::barycentric_coordinates(vList[0],
                                                    vList[1], vList[2], projectedQuery);
  // effectively a clamp, forcing negative bary coords to 0
  int loopIndex = 0;
  double total;
  // could maybe use fewer lines of code passing weight by address
  for (double weight: projected_barys) {
    if (weight < 0) projected_barys[loopIndex] = 0;
    total += projected_barys[loopIndex];
    loopIndex++;
  }
  for (int i = 0; i < 3; i++) {
    baryCoords[i] = projected_barys[i]/total;
  }
}

std::pair<Face_index,Vector_3> selectFaceFromVertex(const Vertex_index &intersectedVertex, const Vector_3 &toIntersection, const Face_index &source_face,
	       const Triangle_mesh &mesh) {
  std::vector<Face_index> candidateFaces; 
  candidateFaces.reserve(10); // usually six, but it can be more  
  
  // get all possible faces that we could walk into: 
  Face_circulator fbegin(mesh.halfedge(intersectedVertex),mesh), done(fbegin); // fbegin is the name of the Face_circulator type object, things 
  									       // in the parentheses are declaration arguments for Face_circulators
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

  Point_3 source_r3 = mesh.point(intersectedVertex);
  
  // below -- we're getting the barycentric coordinates of the vertex by assigning a 1 to the 
  // associated position in an array of otherwise 0s. There is a CGAL function that outputs this, 
  // but it requires me to use a boost graph descriptor wrapper instead of a vertex_index wrapper
  Triangle_mesh::Halfedge_index hf = mesh.halfedge(source_face);
  std::array<double, 3> vert_bary_weights = {0,0,0}; 
  
  int counter = 0; 
  for(Triangle_mesh::Halfedge_index hi : halfedges_around_face(hf, mesh))
  {
    if (source(hi, mesh)==intersectedVertex) {
	vert_bary_weights[counter] = 1;
	break; 
    }
    counter++;  
  }

  Face_location intersectedBary = std::make_pair(source_face, vert_bary_weights); // the face_location of a 
  										  //vertex will be a 1 and two 0 s   

  // now -- check every candidate face to see if the path falls within it after bending. 
  // hypothesis: the path will only fall within the candidate face if that face is the *correct* face, most of the time -- 
  // pathological meshes where two faces are at an extreme angle relative to each other & very thin can break this. 
  // We can check by extending the path out and seeing if its endpoint either intersects the candidate face or falls within
  // that face. If either are true, this is the correct face and we can break the loop.
  Face_index correctFace = source_face;
  bool outOfTriangle;
  Triangle_3 f_triangle; 
  std::array<double, 3> source_bary = intersectedBary.second;
  double closest = 1;
  double near_dist;  
  Face_index closestFace = correctFace;
  Barycentric_coordinates storeCoords;

  Vector_3 bump; 
  
  // we could while loop instead of using break;, but for loop guaranteed terminates
  for (Face_index candidate_face: candidateFaces) {
    std::vector<Vertex_index> faceVertices = getVertexIndices(mesh, candidate_face);
    std::vector<Point_3> fPositions = getVertexPositions(mesh, candidate_face);  
    Face_location candidateLocation = rotateIntoNewFace(mesh, source_face, candidate_face,
		   		    source_r3, source_r3 +  toIntersection);
    f_triangle = Triangle_3(fPositions[0],fPositions[1],fPositions[2]);

    std::array<double, 3> cand_bary = candidateLocation.second;
    
    Point_3 cand_r3 = PMP::construct_point(candidateLocation, mesh); 
 
    near_dist = CGAL::squared_distance(cand_r3,f_triangle); 
    if (near_dist < closest) { 
      closest = near_dist; 
      closestFace = candidate_face;
      storeCoords = cand_bary;
    }

    outOfTriangle = false;
    // if we have points that rotate to be exceptionally close to an edge, we sometimes need 
    // to pick the target that was closest to its candidate face
    for (int i = 0; i < 3; i++) {
      if (cand_bary[i] < 0) outOfTriangle = true;
    } 
     
    if (!outOfTriangle) {
      correctFace = candidate_face;
      bump = normalizer(Vector_3(source_r3, cand_r3)); 
      break;
    }

    if (outOfTriangle) {      
      std::pair<Point_3,std::vector<Vertex_index>> ifIntersect = find_intersection(mesh, source_face, 
		      						    source_r3, cand_r3, faceVertices);
      Point_3 intersection_point = ifIntersect.first; 
      if (intersection_point != Point_3(1000,1000,1000)) { 
	correctFace = candidate_face;
        bump = normalizer(Vector_3(source_r3, cand_r3)); 
	break; 
      }
    } 
  }
  
  // we should always find the right face here, given that we're staying on the mesh.
  if (correctFace == source_face) {
    correctFace = closestFace; 
    std::cout << "No perfectly correct face found, using closest face instead. ";
    std::cout << "\nClosest face: " << closestFace; 
    std::cout << "\nTarget distance to closest face: " << closest << std::endl; 
    
    Face_location finalFLoc = std::make_pair(correctFace,storeCoords);
    Point_3 p1 = PMP::construct_point(finalFLoc,mesh); 
    clampBary(mesh, correctFace, storeCoords);
    Face_location clampedFLoc = std::make_pair(correctFace, storeCoords); 
    Point_3 p2 = PMP::construct_point(clampedFLoc, mesh); 
    Vector_3 v1 = normalizer(Vector_3(source_r3, p1));
    Vector_3 v2 = normalizer(Vector_3(source_r3, p2)); 
    bump = v2; // this should be v2-v1, but for now, want to compare with by angle approach 
  }	

  //std::cout << "Using throughvertex routine, found target face: " << correctFace << std::endl;
  return std::make_pair(correctFace,bump);
}


std::pair<Face_index, Vector_3> throughVertexByAngle(const Vertex_index intersectedVertex, const Vector_3 &toIntersection, const Face_index &source_face,
               const Triangle_mesh &mesh) {
    
    //std::cout << "vertices around vertex " << intersectedVertex << std::endl;
    
    Point_3 vi_r3 = mesh.point(intersectedVertex);
   
    //determine what vertex should be the starting vertex for the angle measurement as 
    //the next vertex from the intersected one (guaranteeing it's in the original face)
    Vertex_index startingVertex;
    
    //if we know what the source face is, this guarantees the first edge
    //we consider will be part of the source face. 
    Halfedge_index hf = mesh.halfedge(source_face);
    for(Halfedge_index hi : halfedges_around_face(hf, mesh))
    {
      if(source(hi, mesh) == intersectedVertex) startingVertex = target(hi,mesh);
    }

    std::vector<Vertex_index> neighborIndices; 
    neighborIndices.reserve(8); 

    Vertex_circulator vbegin(mesh.halfedge(intersectedVertex),mesh), done(vbegin);
    do {
      neighborIndices.push_back(*vbegin++);
    } while(vbegin != done);

    // get vertices as r3 positions starting from intended starting vertex
    std::vector<Point_3> neighborVertices;  
    std::vector<Vertex_index> neighborIndicesCopy; //going to sort this one by starting point first 
    neighborVertices.reserve(neighborIndices.size());
    neighborIndicesCopy.reserve(neighborIndices.size());
    Vertex_index currentVert;
    bool started = false;
    bool atStartPoint = false; 
    bool unfinished = true;
    int ii = -1;  
    while (unfinished) {
      //ii used to iterate through list   
      ii = (ii+1)%neighborIndices.size(); 
      currentVert = neighborIndices[ii];
      atStartPoint = (currentVert == startingVertex);
      if (started && atStartPoint) unfinished=false; 
      else {
        if (atStartPoint) started=true;
        if (started) {
		neighborIndicesCopy.push_back(currentVert); 
		neighborVertices.push_back(mesh.point(currentVert));
        };
      };
    }
    // neighborIndicesCopy stores the vertices in edge traversal order --  
    // coincidentally, therefore, it also stores the ``edges'' as vertex index
    // pairs with the intersected vertex

    // write std::vector of vector objects pointing from the vertex to 
    // each of its neighbors in sequence
    std::vector<Vector_3> edgeVectors;   
    int numNeighbors = neighborVertices.size(); 
    edgeVectors.reserve(numNeighbors); 
    
    for (int i; i < numNeighbors; i++) { 	
	Point_3 v = neighborVertices[i];
        edgeVectors.push_back(normalizer(Vector_3(vi_r3, v)));
    }

    double totalAngle = 0;
    for (int i =0; i < numNeighbors; i++)
        { 
        Vector_3 v1 = edgeVectors[i];
        Vector_3 v2 = edgeVectors[(i+1)%numNeighbors];
        totalAngle += angleBetween(v1,v2);
	} 
    //std::cout << "total angle " << totalAngle << std::endl;

    //now -- pick the vector in a target face that corresponds to having 
    //exactly half of the total angle of the vertex between it and the 
    //vector from the original source to the vertex itself. 
    double angleToTravel = totalAngle/2.0;
    Vector_3 queryHeading = normalizer(-toIntersection); // negative so it matches heading of edge vectors
    double firstAngle = angleBetween(queryHeading,edgeVectors[0]); 
    double angleTraveled = firstAngle;
    int counter = 0; 
    //std::cout << "query heading: " << queryHeading << std::endl;
    //std::cout << "first edge: " << edgeVectors[0] << std::endl;
    //std::cout << "first angle to travel is: " << firstAngle << std::endl; 
    // below should always terminate well  
    while (angleTraveled < angleToTravel) { 
	Vector_3 e1 = edgeVectors[counter]; 
	Vector_3 e2 = edgeVectors[(counter+1)%numNeighbors]; //leaving this here for instances where one face is gigantic and the others are miniscule
        angleTraveled += angleBetween(e1, e2);
        counter ++; // counter will always finish on the e2 index
        if (counter >= numNeighbors) { 
            printf("did not find an angular halfway when picking heading from a vertex, debug!");
	    break;  
	}	
    }
    double angleDifference = angleToTravel - angleTraveled; 
    //std::cout << "after angle traversal, angular distance remaining is " << angleDifference << " and counter is at " << counter << "." << std::endl;
    //std::cout << "angle between last two edges was " << angleBetween(edgeVectors[counter-1],edgeVectors[counter]) << std::endl;
    Vector_3 rotAxisVec = normalizer(CGAL::cross_product(edgeVectors[counter-1], edgeVectors[counter]));
    std::vector<Point_3> rotAxis; 
    rotAxis.reserve(2); 
    rotAxis.push_back(vi_r3);
    rotAxis.push_back(vi_r3 + rotAxisVec); 
    Point_3 rotationTarget = vi_r3 + edgeVectors[counter];  

    Point_3 prospPoint = rotateAboutAxis(rotationTarget, rotAxis, angleDifference);
    Vector_3 prospHeading = normalizer(Vector_3(vi_r3,prospPoint));
    //std::cout << "prospective heading is " << prospHeading << std::endl;
    //std::cout << "angle between prospective heading and first edge: " << angleBetween(edgeVectors[counter-1],prospHeading) << std::endl;
    //std::cout << "and second edge: " << angleBetween(prospHeading,edgeVectors[counter]) << std::endl;
    // now -- need to select target face as face subtended by e1 and e2
    std::vector<Face_index> candidateFaces;
    candidateFaces.reserve(8); // should be six, but managing pathological cases
    Face_circulator fbegin(mesh.halfedge(intersectedVertex),mesh), fdone(fbegin); 
    Face_index cand;
    do {
      cand = *fbegin++; //the next face to add is the face pointed to by the next
                        //iteration of the iterator
      if (cand == source_face) continue; //don't need to consider the face we're coming from
      else {
        candidateFaces.push_back(cand);
      }
    } while(fbegin != fdone);    
    
    Face_index targetFace = source_face;
    std::set<Vertex_index> targetFaceVerts;  
    targetFaceVerts.insert(intersectedVertex);
    targetFaceVerts.insert(neighborIndicesCopy[counter-1]);
    targetFaceVerts.insert(neighborIndicesCopy[counter]);  
    //std::cout << "Target face verts size: " << targetFaceVerts.size() << std::endl;

    for (Face_index face: candidateFaces) { 
        std::vector<Vertex_index> currentFaceVs = getVertexIndices(mesh, face);
        std::set<Vertex_index> currentVsSet; 
	currentVsSet.insert(currentFaceVs[0]);
	currentVsSet.insert(currentFaceVs[1]); 
	currentVsSet.insert(currentFaceVs[2]);
        bool allThere = true;
        for (Vertex_index vert: targetFaceVerts) { 
            auto it = currentVsSet.find(vert);
	    if (it == currentVsSet.end()) 
	        {  
		allThere = false;
		break;
		};
            };
        if (allThere == true)
	    {  
	    targetFace=face;
	    break;
	    }
    } 
    //std::cout << "for source face " << source_face << " found target face " << targetFace << std::endl; 
    //std::cout << "final heading " << prospHeading << std::endl;
    if (source_face == targetFace) throw std::exception(); 
    return std::make_pair(targetFace,prospHeading); 
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

std::pair<Face_index,Vector_3> getTargetFace(std::vector<Vertex_index> intersected, const Vector_3 &toIntersection, const Face_index &source_face, const Triangle_mesh &mesh) {
 //Find the face we're moving into by evaluating the other face attached to the edge made of vertex indices "intersected."
  
  std::string printvar;
  std::cout << "through vertex? ";
  bool throughVertex = (intersected.size() == 1); 
  if (throughVertex) {
    printf("true\n");
    Vertex_index intersectedVertex = intersected[0];
    // usually, below picks from six possible faces. 
    return selectFaceFromVertex(intersectedVertex, toIntersection, source_face, mesh);   
  }
  else printf("false"); 
  Halfedge_index intersected_edge = mesh.halfedge(intersected[0],intersected[1]);
  Face_index provisional_target_face = mesh.face(intersected_edge);
  if (provisional_target_face == source_face) provisional_target_face = mesh.face(mesh.opposite(intersected_edge));
  Vector_3 bump = Vector_3(0,0,0);  
  // ``false'' is a flag to indicate where we're going to go into a face that the rotation will land us outside of
  return std::make_pair(provisional_target_face, bump);
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
  std::pair<Face_index, Vector_3> targetBumpPair;
  Face_index currentTargetFace;
  Vector_3 bump;
  //Vector_3 currentTargetFaceNormal;
  double lengthToSharedElement;
  double rotationAngle;
  double overlap;

  while(intersection){
    //vertexList = getVertexPositions(mesh,currentSourceFace);
    vertexList = getVertexIndices(mesh,currentSourceFace);


    std::pair<Point_3, std::vector<Vertex_index>> intersection_info = find_intersection(mesh, currentSourceFace, source_point,
		                                                                        source_point+current_move, vertexList);
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
    
    targetBumpPair = getTargetFace(intersected_elements, vector_to_intersection, currentSourceFace, mesh); //face we're about to walk into;
    currentTargetFace = targetBumpPair.first;
    bump = targetBumpPair.second;

    source_point = intersection_point;//update source to be the most recent intersection point -- finish walking there 

    Face_location newMoveLocation = rotateIntoNewFace(mesh, currentSourceFace, currentTargetFace, source_point, target);

    std::array<double,3> newBarys = newMoveLocation.second; 

    std::cout << "barycentric coordinates of new position in face " << newMoveLocation.first<< std::endl;
    for (int i = 0; i < 3; i++) std::cout << newBarys[i] << ", ";
    std::cout << std::endl; 
    // if we happen to have rotated into a face without actually being in that face and we know we were supposed to... well, this is difficult. Because
    // then in the error case, for the longest possible shift, we will accidentally cut the shift much, much too short. So just bumping the apparent target 
    // into the face doesn't work. Maybe I can tilt the move vector a little bit to be in the proper face.  
    

    Point_3 rotatedTarget = PMP::construct_point(newMoveLocation,mesh); 

    current_move = Vector_3(source_point, rotatedTarget);//source is now intersection
    // in the very rare instance we're rotating imperfectly into a face -- i.e., the move vector won't 
    // be in it -- we bump the vector into it based on what we learned in face identifier
    if (bump != Vector_3(0,0,0)) { 
      double moveLength = vectorMagnitude(current_move);
      Vector_3 normedMove = normalizer(current_move); 
      current_move = moveLength*(normedMove+bump);
    }

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


