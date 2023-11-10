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


//primary functions
std::pair<std::vector<double>,std::vector<Vector_3>> calcTangentsAndDistances (
		Triangle_mesh mesh, Point_3 source, std::vector<Point_3> targets, std::size_t num_targets) {
  //std::cout << "calculating tangents and distances for source "<< source << std::endl;
  Surface_mesh_shortest_path shortest_paths(mesh);
  AABB_tree tree;
  shortest_paths.build_aabb_tree(tree);
  //convert source point to barycentric coordinates via locate
  const Point_3 source_pt = source;
  Face_location source_loc = shortest_paths.locate<AABB_face_graph_traits>(source_pt,tree);
  shortest_paths.add_source_point(source_loc.first,source_loc.second);
  
  std::vector<double> distances;
  std::vector<Vector_3> tangents;
  std::vector<Point_3> points; 
  for (std::size_t i = 0; i < num_targets; i++) {
    // std::cout << "iteration point " << i+1 << ", at position " << targets[i] <<  "; calling locate" << std::endl;
    Face_location target_loc = shortest_paths.locate<AABB_face_graph_traits>(targets[i],tree);
    //std::cout << "target loc elements: " << target_loc.first << "--" << target_loc.second[0] << " " << target_loc.second[1] << " " << target_loc.second[2] << std::endl;
    shortest_paths.shortest_path_points_to_source_points(target_loc.first, target_loc.second, std::back_inserter(points));
    //std::cout << "calculating shortest distance" << std::endl;
    distances.push_back(std::get<0>( shortest_paths.shortest_distance_to_source_points(target_loc.first, target_loc.second)));
    //std::cout << "distance calculated." << std::endl;
    //path goes from target to source -- so if we want to know path tangent at 
    //source for force calculation, we must use the *end* of points[]
    tangents.push_back(Vector_3(points[points.size()-2],points[points.size()-1]));
    points.clear();
  }
 
  return std::make_pair(distances,tangents);
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

Vector_3 force_on_source (Triangle_mesh mesh, Point_3 source, std::vector<Point_3> targets, std::size_t num_targets) {
  //create list of distances and path tangents between the source particle and the targets
  std::pair<std::vector<double>, std::vector<Vector_3>> distancesAndTangents = calcTangentsAndDistances(mesh, source, targets, num_targets);   
  std::vector<double> distances = distancesAndTangents.first;
  std::vector<Vector_3> tangents = distancesAndTangents.second;
  //we could either do this loop within forceFunction or here -- shouldn't be hugely different
  
  //L-J parameters
  double epsilon = 1;
  double sigma = 1.1;
  Vector_3 force= Vector_3(0,0,0); //initialize to zero to avoid redefinition --
                                   //also handles case of no neighbors
				   
  for (std::size_t i = 0; i < distances.size(); i++) {
    force+= LJForce(distances[i],tangents[i], epsilon, sigma);    
  } 
  std::cout << "Calculated force magnitude is " << vectorMagnitude(force) << std::endl; 
  return force;
}

