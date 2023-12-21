//Author: Toler H. Webb
//Description: CGAL code to calculate the force between two particles at given locations. 
//             Particle shape information is stored in the force function. 

//includes:
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include "utils.h"

//typedefs:
typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     KernelWithSqrt;
typedef Kernel::Point_3                                                 Point_3;
typedef Kernel::Vector_3                                                Vector_3;
typedef Kernel::Ray_3                                                   Ray_3;
typedef CGAL::Surface_mesh<Point_3>                                     Triangle_mesh;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh>  Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits>                        Surface_mesh_shortest_path;
typedef boost::graph_traits<Triangle_mesh>                              Graph_traits;
typedef Graph_traits::vertex_iterator                                   vertex_iterator;
typedef Graph_traits::face_iterator                                     face_iterator;
typedef typename Surface_mesh_shortest_path::Barycentric_coordinates    Barycentric_coordinates;
typedef typename Surface_mesh_shortest_path::Face_location              Face_location;
typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh>         AABB_face_graph_primitive;
typedef CGAL::AABB_traits<Kernel, AABB_face_graph_primitive>            AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                         AABB_tree;


int time_big_sequence_tree ( const Triangle_mesh &mesh, const std::vector<Point_3> &source) {
  std::chrono::duration<long, std::milli> seq_time;
  Surface_mesh_shortest_path shortest_paths(mesh);
  AABB_tree tree;
  shortest_paths.build_aabb_tree(tree);
  Face_location source_loc;

  for (Point_3 sp: source) {
    source_loc = shortest_paths.locate<AABB_face_graph_traits>(sp,tree);
    shortest_paths.add_source_point(source_loc.first,source_loc.second);
  }  	  

  auto start = std::chrono::high_resolution_clock::now();
  shortest_paths.build_sequence_tree();
  auto end = std::chrono::high_resolution_clock::now();
  seq_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "time to make tree with " << source.size() << " source particles: " << seq_time.count() << "ms" << std::endl;

  return seq_time.count();
}

//primary functions
std::pair<std::vector<double>,std::vector<Vector_3>> calcTangentsAndDistances (
		const Triangle_mesh &mesh, const Point_3 &source, const std::vector<Point_3> &targets, const std::size_t &num_targets) {
 
  std::chrono::duration<long, std::milli> f_time;
  auto overstart = std::chrono::high_resolution_clock::now();

  auto start = std::chrono::high_resolution_clock::now();
  Surface_mesh_shortest_path shortest_paths(mesh);
  AABB_tree tree;
  shortest_paths.build_aabb_tree(tree);
  auto end = std::chrono::high_resolution_clock::now();
  f_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "instantiate shortest paths & build aabb tree: " << f_time.count() << "ms" <<  std::endl;
  
  //convert source point to barycentric coordinates via locate
  start = std::chrono::high_resolution_clock::now();
  const Point_3 source_pt = source;
  Face_location source_loc = shortest_paths.locate<AABB_face_graph_traits>(source_pt,tree);
  end = std::chrono::high_resolution_clock::now();
  f_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "initial locate: " << f_time.count() << "ms" <<  std::endl;
  
  start = std::chrono::high_resolution_clock::now();
  shortest_paths.add_source_point(source_loc.first,source_loc.second);
  end = std::chrono::high_resolution_clock::now();
  f_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "time to add source point: " << f_time.count() << "ms" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  shortest_paths.build_sequence_tree();
  end = std::chrono::high_resolution_clock::now();
  f_time =  std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "time to build sequence tree: " << f_time.count() << "ms" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  std::vector<double> distances;
  std::vector<Vector_3> tangents;
  std::vector<Point_3> points; 
  //"distances" and "tangents" will both store data from a path between a 
  // pair of particles, so we can directly reserve # of pairs
  distances.reserve(num_targets);
  tangents.reserve(num_targets);
  end = std::chrono::high_resolution_clock::now();
  f_time =  std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "instantiate storage vectors: " << f_time.count() << "ms" <<  std::endl;


  start = std::chrono::high_resolution_clock::now();
  for (std::size_t i = 0; i < num_targets; i++) {
    Face_location target_loc = shortest_paths.locate<AABB_face_graph_traits>(targets[i],tree);
    
    //start = std::chrono::high_resolution_clock::now();
    shortest_paths.shortest_path_points_to_source_points(target_loc.first, target_loc.second, std::back_inserter(points));
    //end = std::chrono::high_resolution_clock::now();
    //f_time =  std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    //std::cout << "time to compute path: " << f_time.count() << std::endl;

    distances.push_back(std::get<0>( shortest_paths.shortest_distance_to_source_points(target_loc.first, target_loc.second)));
    //path goes from target to source -- so if we want to know path tangent at 
    //source for force calculation, we must use the *end* of points[]
    tangents.push_back(Vector_3(points[points.size()-2],points[points.size()-1]));
    points.clear();
  }
  end = std::chrono::high_resolution_clock::now();
  f_time =  std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "for loop: " << f_time.count() << "ms" <<  std::endl;

  auto overend = std::chrono::high_resolution_clock::now();
  f_time = std::chrono::duration_cast<std::chrono::milliseconds>(overend-overstart);
  std::cout << "in-function total time: " << f_time.count() << "ms" << std::endl;

  start = std::chrono::high_resolution_clock::now(); 
  std::pair<std::vector<double>,std::vector<Vector_3>> pairforreturn = std::make_pair(distances,tangents);
  end = std::chrono::high_resolution_clock::now();
  f_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "make_pair at the end: " << f_time.count() << "ms" << std::endl;


  return pairforreturn;
}


Vector_3 LJForce (float dist, Vector_3 tangent, double epsilon, double sigma) { 
  //once again, a=1 and c=3
  Vector_3 normalizedTangent = normalizer(tangent);

  //lennard-jones potential is 
  //4 epsilon((sigma/r)^12 - (sigma/r)^6) 
  //associated force is 
  //4 epsilon (12 sigma^12 / r^13 - 6 sigma ^6 / r^7)nabla(r)
  //where r is the distance function on the surface 
  
  Vector_3 force = 4*epsilon*(12*pow(sigma,12)/pow(dist,13)-6*pow(sigma,6)/pow(dist,7))*normalizedTangent;
  
  return force;
}


Vector_3 simpleRepulsion(float dist, Vector_3 tangent, double sigma) {
  //simple repulsive force from a Gaussian potential for easy checking.  
  
  //for me, writing below as literally negative gradient of Gaussian potential. 
  //maximum value is (2^(3/2))/sigma)e^(-2). 

  Vector_3 normalizedTangent = normalizer(tangent); 
  Vector_3 force = -(-2*dist/pow(sigma,2))*exp(-pow(dist,2)/pow(sigma,2))*normalizedTangent;

  return force;

}


Vector_3 inversePowerLaw(float dist, Vector_3 tangent, double sigma) {
 //fill in later
 return Vector_3(0,0,0);
}

Vector_3 force_on_source (const Triangle_mesh &mesh, const Point_3 &source, const std::vector<Point_3> &targets, const std::size_t &num_targets) {
  //create list of distances and path tangents between the source particle and the targets
  std::chrono::duration<long, std::milli> f_time;
  auto start = std::chrono::high_resolution_clock::now();
  std::pair<std::vector<double>, std::vector<Vector_3>> distancesAndTangents = calcTangentsAndDistances(mesh, source, targets, num_targets);   
  auto end = std::chrono::high_resolution_clock::now();
  f_time =  std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "time to calc dist & tang:  " << f_time.count() << std::endl;
  std::cout << "" << std::endl;

  std::vector<double> distances = distancesAndTangents.first;
  std::vector<Vector_3> tangents = distancesAndTangents.second;
  //we could either do this loop within forceFunction or here -- shouldn't be hugely different
  
  //L-J parameters
  double epsilon = 1;
  double sigma = 1;
  Vector_3 force= Vector_3(0,0,0); //initialize to zero to avoid redefinition --
                                   //also handles case of no neighbors
				   
  for (std::size_t i = 0; i < distances.size(); i++) {
    force+= simpleRepulsion(distances[i],tangents[i], sigma);
  } 
  //std::cout << "Calculated force magnitude is " << vectorMagnitude(force) << std::endl; 
  return force;
}

/**
 * Protective function to ensure that we don't add duplicate vertices to the mesh. 
 * This would be really costly if submesh got really large (as every time we added a new vertex,
 * we would need to look through *all* the vertices) but, in principle, it shouldn't get big enough
 * to matter.  
 */
Vertex_index unique_vertex_to_mesh(Triangle_mesh& mesh, const Point_3 new_vertex) {
    // Check if the vertex already exists in the mesh
    Vertex_index vi = Triangle_mesh::null_vertex();
    for (Vertex_index v : mesh.vertices()) {
        if (mesh.point(v) == new_vertex) {
	    vi = v;
            break;
        }
    }
    // If the vertex wasn't found in the above loop, create a new vertex and add it to the mesh
    if (vi == Triangle_mesh::null_vertex()) {
      vi = mesh.add_vertex(new_vertex); 
    } 

    // Return the index of the intended vertex so we can add faces
    return vi; 
}  

Triangle_mesh create_submesh_from_visited(std::set<Face_index> visited, const Triangle_mesh& mesh) { 
  Triangle_mesh submesh; 
  
  Vertex_index u; 
  Vertex_index v;
  Vertex_index w; 
  std::vector<Point_3> face_vertices; 
  Face_index added_face; 
  // because we use a face_index *set*, face indices in visited are guaranteed to be unique.
  // Thus, we shouldn't need to check for duplicate face additions. However, there's no 
  // such guarantee for the vertices.  
  for (Face_index face: visited) { 
    face_vertices = getVertexPositions(mesh, face);
    Vertex_index u = unique_vertex_to_mesh(submesh, face_vertices[0]);
    Vertex_index v = unique_vertex_to_mesh(submesh, face_vertices[1]);
    Vertex_index w = unique_vertex_to_mesh(submesh, face_vertices[2]);
    added_face = submesh.add_face(u,v,w); 
    // below to make sure we haven't borked the mesh to start with
    if(added_face == CGAL::Surface_mesh<Kernel>::null_face())
    {
      std::cerr<<"Orientation error in vertices"<<std::endl;
      added_face = submesh.add_face(u,w,v);
      assert(added_face != Triangle_mesh::null_face());
    }
  }

  return submesh; 
}

/**
 * Code to generate the submesh of minimum size that includes source, targets, and all intermediary faces.
 * 
 */
Triangle_mesh build_minimum_submesh(const Face_location& source, const std::vector<Face_location>& targets,
                                    const double& cutoff_dist, const Triangle_mesh& mesh) {
  // Face_location objects are std::pairs of face indices and barycentric coordinates
  auto start = std::chrono::high_resolution_clock::now();
	
  Triangle_mesh submesh; 
  Point_3 source_r3 = PMP::construct_point(source,mesh); 

  Face_index source_face = source.first;
   
  std::set<Face_index> goal_faces; 
  
  std::set<Face_index> visited; 
  visited.insert(source_face);

  // first: add to the list of goal faces all the faces which contain 
  // a target but are not the source face
  for (Face_location target: targets) { 
    if (target.first != source_face) goal_faces.insert(target.first); 
  }
  if (goal_faces.empty()) return create_submesh_from_visited(visited,mesh); // early return if all our targets are in one face

  std::vector<Face_index> exploration_stack; 
  std::vector<Point_3> face_verts; 
  face_verts.reserve(3);
  Face_index neighboring_face;
  
  Halfedge_index hf = mesh.halfedge(source_face);

  for(Halfedge_index hi : halfedges_around_face(hf, mesh)){
    // Access the neighboring face through the opposite halfedge of the current halfedge
    neighboring_face = mesh.face(mesh.opposite(hi));
    face_verts = getVertexPositions(mesh, neighboring_face);
    visited.insert(neighboring_face); 
    exploration_stack.push_back(neighboring_face);
    
    if(goal_faces.count(neighboring_face)) {
      goal_faces.erase(neighboring_face);
      // since we know there can't be holes in a mesh with only faces that are directly connected 
      // to the source face: 
      if (goal_faces.empty()) {
        return create_submesh_from_visited(visited, mesh); 
      } 
    } 
    // push to queue to investigate neighbors so long as we still have targets to find
  }

  Face_index current_face;
  // now: breadth-first search, checking every connected face and trying to make a complete patch
  while (!exploration_stack.empty() && !goal_faces.empty()) { 
    // get the current face we're "standing" in
    current_face = exploration_stack.back(); 
    exploration_stack.pop_back(); 

    // redefining hf since it will be serving the same purpose
    hf = mesh.halfedge(current_face); 

    for (Halfedge_index hi : halfedges_around_face(hf, mesh)) {
      neighboring_face = mesh.face(mesh.opposite(hi));
      // if we've seen this face before, continue so we don't create a loop
      if(visited.count(neighboring_face)) continue;
      // otherwise, say we've been here before
      
      face_verts = getVertexPositions(mesh,neighboring_face);
      // if all our vertices are outside the cutoff radius, we're out of bounds
      // and shouldn't add this face to the exploration stack
      bool allOut = true;
      for (Point_3 v : face_verts) {
         if (vectorMagnitude(Vector_3(source_r3,v)) < cutoff_dist) {
           allOut = false;
	   break;
	   }
      }
      if (allOut) continue;
       
      // the two statements above guarantee closure, although closed surfaces 
      // with small internal "radii"/large cutoff radii can lead to delayed closure/
      // submeshes that are too large becuase they will require the first condition

      visited.insert(neighboring_face); 
      
      // for now, this is just a debug statement -- make sure that we've
      // visited all goal faces later... but we should be able to continue
      // to reduce the submesh size, no? I'm worried about holes, but that's it. 
      if(goal_faces.count(neighboring_face)) {
        goal_faces.erase(neighboring_face); 
	continue;
      }
      
      //finally, if we absolutely need to keep looking down this path (it is not a goal,
      //it is not past an edge, we haven't seen it before), we do. 
      exploration_stack.push_back(neighboring_face); 
    }    

  }
  if (!goal_faces.empty()) {
    std::cout << "All goal faces not found, debug!" << std::endl;	
  }

  auto end = std::chrono::high_resolution_clock::now(); 
  auto  traversal_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start); 
  std::cout << "time to execute mesh traversal for submesh generation: " << traversal_time.count() << " microseconds" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  submesh = create_submesh_from_visited(visited,mesh); 
  end = std::chrono::high_resolution_clock::now();
  auto submesh_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
  std::cout << "time to create submesh from visited: "  << submesh_time.count() << " microseconds" << std::endl;
  return create_submesh_from_visited(visited,mesh); 
}

