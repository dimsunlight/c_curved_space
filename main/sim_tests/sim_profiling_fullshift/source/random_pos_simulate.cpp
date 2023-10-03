/* Author: Toler H. Webb
 * Description: CGAL-based routine for simulating the interaction of a 
 * small number of particles.  
 */

//includes
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
#include <CGAL/Point_set_3.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Point_set_3/IO/XYZ.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <math.h>
#include "shift.h"
#include "get_force.h"
#include <chrono>
#include <ratio>
#include <thread>
#include <random>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     KernelWithSqrt;
typedef Kernel::Point_3                                                 Point_3;
typedef Kernel::Vector_3                                                Vector_3;
typedef CGAL::Point_set_3<Point_3>                                      Point_set_3;
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


float random_in_range(float range_min, float range_max)
{
  //uses mersenne twister to get a random real #
  std::random_device                   rand_dev;
  std::mt19937                         generator(rand_dev());
  std::uniform_real_distribution<float>  distr(range_min, range_max);

  return distr(generator);
}

float torus_weight(float a, float c, float u, float v) {
  //get torus weight for sampling against [0,1] distribution --
  //so we divide by radius to compare to ``unit'' torus
  return (c+a*cos(v))/(c+a);
}

//float sinplane_weight --- look up what differential geometric aspects are needed for rejections ampling

Point_3 torus_sample(float a, float c) {
  //rejection method sample on a torus; that is, reject anything without the appropriate weight
  float U;
  float V;
  float W;
  float weight;
  Point_3 returnPoint = Point_3(NULL,NULL,NULL);
  while (returnPoint == Point_3(NULL,NULL,NULL)) {
    U = random_in_range(0.0, 2*M_PI);
    V = random_in_range(0.0, 2*M_PI);
    W = random_in_range(0.0, 1.0);
    weight = torus_weight(a,c,U,V);
    if (W <= weight) {
      returnPoint = Point_3((c+a*cos(V))*cos(U), (c+a*cos(V))*sin(U), a*sin(V));
      return returnPoint;
    }
  }
  //this should never happen, but just to avoid false error messages
  return Point_3(NULL,NULL,NULL);//filler point to check if we have sampled at all
}

//Point_3 sinplane_sample() {} ; once weight is defined, write as above

std::vector<Point_3> n_torus_sample_points(std::size_t n, float a, float c) {
  //I use std::vectors of point_3 objects for positions, so we loop
  //to create a similar vector for n random positions
  std::vector<Point_3> sample_points = {};
  for (std::size_t i=0; i < n; i++) {
    sample_points.push_back(torus_sample(a,c));
  }
  return sample_points;  
}


std::vector<Point_3> create_particles_from_xyz(std::string locations_file) {
  Point_set_3 locations;
  CGAL::IO::read_XYZ(locations_file,locations);
  std::vector<Point_3> outputPoints; //I should probably be using
                                    //Point_set_3 objects throughout my code; cleanup for later
  for (Point_set_3::const_iterator it = locations.begin(); it != locations.end(); ++ it) {
    outputPoints.push_back(locations.point(*it));
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


float averageTime(std::vector<std::chrono::milliseconds> times,std::size_t timesteps) {
  int sum = 0;
  float n=0.0;
  for (std::chrono::milliseconds element: times) {
    sum += element.count();
    n+= 1.0/timesteps; //scaled as the length of 1 millisecond relative to a second 
            //(n+=1 would be adding a second to the division instead of a millisecond)
  }

  return sum/n;
}

int main (int argc, char* argv[]) {

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("sims_project/torusrb20.off");
  Triangle_mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  std::size_t timesteps = 10;
  double      stepsize  =.001;
  double  neighbor_cutoff = 8;
  std::size_t n_simulations = 50;

  std::size_t num_particles;
  std::vector<Point_3> particle_locations;
  std::vector<Point_3> location_buffer;
  std::vector<std::chrono::milliseconds> timestep_costs;
  std::vector<std::pair<Point_3,std::vector<Point_3>>> particles_with_neighbors;
  std::pair<Point_3,std::vector<Point_3>> particle_and_neighbors;
  std::chrono::steady_clock::time_point sim_start;
  std::chrono::steady_clock::time_point sim_end;
  std::chrono::steady_clock::time_point timestep_start;
  std::chrono::steady_clock::time_point timestep_end;
  Vector_3 f_on_p;

  const std::string times_filename = "torus_times.txt";
  std::ofstream times_file(times_filename);
  std::cout << "Making simulations, timestep size " << stepsize << " for " << timesteps << "timesteps." <<  std::endl;
  
  for (std::size_t iter = 0; iter < n_simulations; iter++) {
    std::cout << "run number: " << iter << std::endl;
    num_particles = 2+iter*4;  
    particle_locations = n_torus_sample_points(num_particles, 1, 3);

    //std::cout << "printing locations: " << std::endl;
    //for (Point_3 location: particle_locations) std::cout << "(" <<location << ")"  << std::endl;

    particles_with_neighbors = get_neighbors(particle_locations,neighbor_cutoff);
  
    sim_start = std::chrono::steady_clock::now();
    //main simulation loop
    for (std::size_t j = 0; j < timesteps; j++) {
      
      timestep_start = std::chrono::steady_clock::now();
      particles_with_neighbors = get_neighbors(particle_locations,neighbor_cutoff);
    
      //find forces and do shift
      for (std::size_t i = 0; i < particle_locations.size();i++) {
        particle_and_neighbors = particles_with_neighbors[i];
        f_on_p = force_on_source(mesh,particle_and_neighbors.first,particle_and_neighbors.second,
		    particle_and_neighbors.second.size()); 
        location_buffer.push_back(shift(mesh, particle_and_neighbors.first, stepsize*f_on_p));
      }
      //location buffer housekeeping
      particle_locations.clear();
      for (Point_3 location: location_buffer) {
        particle_locations.push_back(location);
      }
      location_buffer.clear();
      timestep_end = std::chrono::steady_clock::now();
      // std::cout << "Timestep cost: " << std::chrono::duration_cast<std::chrono::milliseconds>(timestep_end-timestep_start).count() << "ms" << std::endl; 
      timestep_costs.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(timestep_end-timestep_start));  
    }
       
    sim_end = std::chrono::steady_clock::now();

    std::cout << "Time taken for simulation with "<< num_particles << " particles: " << std::chrono::duration_cast<std::chrono::milliseconds>(sim_end - sim_start).count() << "ms" << std::endl;
    std::cout << "Average time taken per timestep: " << averageTime(timestep_costs,timesteps)/timesteps << "ms" <<  std::endl;
    times_file << num_particles << ", " << averageTime(timestep_costs,timesteps)/timesteps << std::endl;
    timestep_costs.clear();
  }

  return 0;
}





