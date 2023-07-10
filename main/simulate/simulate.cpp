/* Author: Toler H. Webb
 * Description: CGAL-based routine for simulating the interaction of a 
 * single particle pair. 
 */

//includes
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
#include "shift.h"
#include "get_force.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     KernelWithSqrt;
typedef Kernel::Point_3                                                 Point_3;
typedef Kernel::Vector_3                                                Vector_3;
typedef Kernel::Point_set_3                                             Point_set_3;
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

//utility functions

float distBetween(Point_3 point1,Point_3 point2) {
  //I could (should) use geodesic distance for this, but euclidean distance
  //is
  //  1) easier to debug
  //  2) faster for now
  //so for prototyping, I'm going to stick with this. 

  //vector between p1 and p2, using v because I have to use accessors 
  Vector_3 v = Vector_3(point1,point2);
  float slength = v.x()*v.x() + v.y()*v.y()+v.z()*v.z();
  
  return sqrt(slength);
}

std::vector<Point_3> create_particles_from_xyz(locations_file) {
  Point_set_3 locations = CGAL::IO::read_xyz(locations_file);
  std::vector<Point_3> outputPoints; //I should probably be using
                                    //Point_set_3 objects throughout my code; cleanup for later
  for (Point_3 i: locations) {
    outputPoints.push_back(i);
  }
  return outputPoints;
}	

std::vector<std::pair<Point_3,std::vector<Point_3>>> get_neighbors(std::vector<Point_3> positions, float cutoff){
  //loop through all positions in get_neighbors
  //for each particle 
  //store particle as source (first half of pair)
  //initialize list of neighbors
  //for every particle within (euclidean) cutoff of current source:
  //	add particle location to list of neighbors
  //return list of (particle, neighbors) pairs.
  Point_3 source;
  std::vector<Point_3> neighbors;
  std::vector<std::pair<Point_3,std::vector<Point_3>>> particles_neighbors_pairs; 
  for (Point_3 i: positions) {
    source = i;
    //look for list comprehension method to use later, rather than loop back through
    for (Point_3 j: positions) {
      if (j != i) {
        if (distBetween(i,j) < cutoff) {
	  neighbors.push_back(j);
	}	
      }
    }
    particles_neighbors_pairs.push_back(std::make_pair(source,neighbors));
    neighbors.clear();
  } 
  
  return particles_neighbors_pairs;
}

int main (int argc, char* argv[]) {

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("sims_project/torusrb20.off");
  Triangle_mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  const std::string loc_filename = (argc > 2) ? argv[2] : CGAL::data_file_path("sims_project/default_pos.xyz");
  std::vector<Point_3> particle_locations = create_particles_from_xyz(loc_filename);

  //have defaults for both loading functionalities below
  particle_locations = create_particles_from_xyz(input_xyz); 
  Point_3 location_buffer[particle_locations.size()];

  std::vector<std::pair<Point_3,std::vector<Point_3>>> particles_with_neighbors = get_neighbors(particle_locations);
  
  //can make this multi-step simply by enclosing into timestep loop
  std::pair<Point_3,std::vector<Point_3>> particle_and_neighbors;
  Vector_3 f_on_p;
  for (std::size_t i = 0; i < location_buffer.size();i++) {
    particle_and_neighbors = particles_with_neighbors[i];
    f_on_p = force_on_source(mesh,particle_and_neighbors.first,particle_and_neighbors.second,
		    particle_and_neighbors.second.size()); 
    locationBuffer[i] = shift(mesh, particle_and_neighbors.first,f_on_p);
  }
  //here, we could store location_buffers in a std::vector to use them
  //for animation later
  particle_locations = location_buffer;
  

  
  return 0;
}





