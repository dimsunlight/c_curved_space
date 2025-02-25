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
#include "utils.h" 
#include <chrono>
#include <ratio>
#include <thread>
#include <random>
#include <omp.h>

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
typedef std::mt19937                                                    MTgenerator; //mersenne-twister generator
typedef std::uniform_real_distribution<float>                           distribution;

//unique utility functions

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

MTgenerator make_generator(int seed)
{
  //create random number generator using seed parameter. 

  std::random_device                  rand_dev;
  MTgenerator                         generator(rand_dev());
  generator.seed(seed); //seeded for reproducibility -- make sure we use the same sequence every time :)  

  return generator;
}

distribution distribution_of_range(float range_min, float range_max) {
  distribution distro(range_min, range_max);
  return distro; 
}

float sinplane_weight(float x, float y) {
  return 1/sqrt(1+cos(x)*cos(x)/4.0 + cos(y)*cos(y)/4.0);
}

float torus_weight(float a, float c, float u, float v) {
  //get torus weight for sampling against [0,1] distribution --
  //so we divide by radius to compare to ``unit'' torus
  return (c+a*cos(v))/(c+a);
}

Point_3 sinplane_sample(MTgenerator& gen) {
  //pass by address above to ensure that we sequentially sample 
  //from the generator without using the same starting point.
  distribution SPdistro = distribution_of_range(-1.0,1.0);
  distribution weightDistro = distribution_of_range(0.0,1.0); 
  float X;
  float Y;
  float W;
  float weight;
  Point_3 returnPoint = Point_3(NULL,NULL,NULL);
  while (returnPoint == Point_3(NULL,NULL,NULL)) {
    //goal -- reference gen by address so each call to gen advances the 
    //rng forward once
    X = SPdistro(gen);
    Y = SPdistro(gen);
    W = weightDistro(gen);
    weight = sinplane_weight(X,Y);
    if (W <= weight) {
      returnPoint = Point_3(X,Y,sin(X));
      return returnPoint;
    }
  }
  return returnPoint;
}

Point_3 torus_sample(float a, float c, MTgenerator& gen) {
  //rejection method sample on a torus; that is, reject anything above the corresponding weight. 
  //we pass a mersenne-twister generator by reference 
  distribution angularDistro = distribution_of_range(0.0, 2*M_PI);
  distribution weightDistro = distribution_of_range(0.0,1.0); 
  float U;
  float V;
  float W;
  float weight;
  Point_3 returnPoint = Point_3(NULL,NULL,NULL);
  while (returnPoint == Point_3(NULL,NULL,NULL)) {
    U = angularDistro(gen);
    V = angularDistro(gen);
    W = weightDistro(gen);
    weight = torus_weight(a,c,U,V);
    if (W <= weight) {
      returnPoint = Point_3((c+a*cos(V))*cos(U), (c+a*cos(V))*sin(U), a*sin(V));
      return returnPoint;
    }
  }
  //this should never happen, but just to avoid false error messages
  return Point_3(NULL,NULL,NULL);//filler point to check if we have sampled at all
}

std::vector<Point_3> n_torus_sample_points(std::size_t n, float a, float c, int seed=1997) {
  //I use std::vectors of point_3 objects for positions, so we loop
  //to create a similar vector for n random positions
  std::vector<Point_3> sample_points = {};
  sample_points.reserve(n);
  MTgenerator gen = make_generator(seed); 
  
  for (std::size_t i=0; i < n; i++) {
    sample_points.push_back(torus_sample(a,c, gen));
  }
  return sample_points;  
}

std::vector<Point_3> n_sinplane_sample_points(std::size_t n, int seed = 1997) {
  std::vector<Point_3> sample_points = {};
  sample_points.reserve(n);
  MTgenerator gen = make_generator(seed); 

  for (std::size_t i=0; i < n; i++) {
    sample_points.push_back(sinplane_sample(gen));
  }
  return sample_points;
}

std::vector<Point_3> create_particles_from_xyz(std::string locations_file) {
  Point_set_3 locations;
  CGAL::IO::read_XYZ(locations_file,locations);
  std::vector<Point_3> outputPoints = {}; //I should probably be using
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


std::chrono::milliseconds averageTime(std::vector<std::chrono::milliseconds> times) {
  std::chrono::milliseconds sum;
  int n=0;
  for (std::chrono::milliseconds element: times) {
    sum += element;
    n+=1;
  }

  return sum/n;
}

int main (int argc, char* argv[]) {

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("torusrb20.off");
  Triangle_mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  const std::string loc_filename = (argc > 2) ? argv[2] : "FILLER.xyz";
  
  std::vector<Point_3> particle_locations; 
  if (loc_filename == "FILLER.xyz") {
    particle_locations = n_torus_sample_points(100, 1, 3); 
  } else { 
    particle_locations = create_particles_from_xyz(loc_filename);
  }
  //have defaults for both loading functionalities below
  std::cout << "Original locations:" << std::endl;
  for (Point_3 location: particle_locations) std::cout << location << ", ";
  std::cout << std::endl;
 
  //simulation time parameters -- too-large step sizes can break shift!
  std::size_t timesteps = 1000; //hard defining this now rather than input-defining to avoid extra debugging 
  double      stepsize = .001;   //smaller stepsize for lj (.00001) -- larger for gaussian repulsion

  //define location buffer to ensure simultaneous position update, define neighbor lists, initalize loop variables
  std::vector<Point_3> location_buffer;
  location_buffer.reserve(particle_locations.size());
  std::cout << "Running simulation with " << particle_locations.size() << " particles and " << timesteps << " timesteps." << std::endl;
  std::cout << "Timestep size " << stepsize << std::endl;
  double neighbor_cutoff = 3;
  std::vector<std::pair<Point_3,std::vector<Point_3>>> particles_with_neighbors = get_neighbors(particle_locations,neighbor_cutoff);
  std::pair<Point_3,std::vector<Point_3>> particle_and_neighbors;
  Vector_3 f_on_p;
  
  //write trajectory data to file
  const std::string trajectory_filename = "torus_GR_1000pts_sim.txt";
  std::ofstream trajectory_file(trajectory_filename);
  if (trajectory_file.is_open()){
    trajectory_file << particle_locations.size();
    trajectory_file << "\n";
    for (Point_3 location: particle_locations) {
      trajectory_file << location;
      trajectory_file << "\n"; 
    }
  }


  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
  
  std::cout << "Time taken for setup: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;
  
  std::vector<std::chrono::milliseconds> forceTimes = {};
  std::vector<std::chrono::milliseconds> shiftTimes = {};
  std::chrono::steady_clock::time_point sim_start = std::chrono::steady_clock::now();
  //main simulation loop
  for (std::size_t j = 0; j < timesteps; j++) {
    std::cout << "timestep " << j << " first ten locations: " << std::endl;
    for (int i = 0; i < 10; i++) std::cout << particle_locations[i] << std::endl;
    particles_with_neighbors = get_neighbors(particle_locations,neighbor_cutoff);
    
    //find forces and do shift -- this is the step that most needs to be parallelized,
    //so we can build many sequence trees simultaneously 
    #pragma omp parallel for 
    for (std::size_t i = 0; i < particle_locations.size();i++) {
      particle_and_neighbors = particles_with_neighbors[i];
      //std::cout << "current particle: " <<particle_and_neighbors.first << std::endl;
      auto f_start = std::chrono::high_resolution_clock::now();
      f_on_p = force_on_source(mesh,particle_and_neighbors.first,particle_and_neighbors.second,
		    particle_and_neighbors.second.size()); 
      auto f_end = std::chrono::high_resolution_clock::now();
      auto forcetime = std::chrono::duration_cast<std::chrono::milliseconds>(f_end - f_start);; 
      forceTimes.push_back(forcetime);


      auto s_start = std::chrono::high_resolution_clock::now(); 
      //std::cout << "The force on the particle at " << particle_and_neighbors.first << " is " << f_on_p << std::endl; 
      location_buffer.push_back(shift(mesh, particle_and_neighbors.first, stepsize*f_on_p));
      auto s_end   = std::chrono::high_resolution_clock::now(); 
      std::chrono::duration<long, std::milli> shifttime =  std::chrono::duration_cast<std::chrono::milliseconds>(s_end-s_start);
      shiftTimes.push_back(shifttime);
    }
    //location buffer housekeeping
    particle_locations.clear();
    
    for (Point_3 location: location_buffer) {
      particle_locations.push_back(location);
      //write every ten steps
      if (j%10 == 0) { 
      	trajectory_file << location;
      	trajectory_file << "\n";
      }	
    }
    location_buffer.clear();
    trajectory_file << "\n";
  }
  trajectory_file.close();
  for (Point_3 location: particle_locations) std::cout << location << std::endl;
  std::chrono::steady_clock::time_point sim_end = std::chrono::steady_clock::now();
  std::cout << "Time taken for simulation: " << std::chrono::duration_cast<std::chrono::milliseconds>(sim_end - sim_start).count() << "ms" << std::endl;
  std::cout << "Time for force calculation: " << averageTime(forceTimes).count() << "ms" << std::endl;
  std::cout << "Time for shift calculation: " << averageTime(shiftTimes).count() << "ms" << std::endl;


  return 0;
}





