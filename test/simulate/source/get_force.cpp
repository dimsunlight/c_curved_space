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

