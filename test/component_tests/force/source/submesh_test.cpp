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
#include "get_force.h"
#include "utils.h"
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
typedef std::mt19937                                                    MTgenerator; //mersenne twister generator
typedef std::uniform_real_distribution<float>                           distribution;

namespace PMP = CGAL::Polygon_mesh_processing; 


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
  float sp_ext = 5.0;
  distribution SPdistro = distribution_of_range(-sp_ext,sp_ext);
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

int main(int argc, char* argv[]) {
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  const std::string read_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("torus_isotropic.off");
  Triangle_mesh mesh;
  
  if(!CGAL::IO::read_polygon_mesh(read_filename,mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // basic test of submeshing, agnostic of other dependencies, force calculation itself, 
  std::cout <<  "Testing submeshing routine with mesh file " << read_filename << std::endl; 
  Point_3 point1;
  std::vector<Point_3> point2; 
  std::chrono::duration<long, std::milli> f_time;
  std::size_t targetcount;
  
  int npoints = 100000;
  
  std::vector<Point_3> particle_xyz_locations = n_torus_sample_points(npoints, 1, 3);
  std::vector<Face_location> target_locs; 
  double cutoff_rad = .5;
  Point_3 source_pt = particle_xyz_locations[0];
  Face_location source_loc = PMP::locate(source_pt, mesh);
  printf("Source pt: "); 
  std::cout << source_pt << std::endl;

  printf("target points: ");
  for (Point_3 pt: particle_xyz_locations) {
    if (pt != source_pt) {
      if (vectorMagnitude(Vector_3(source_pt, pt)) < cutoff_rad) {
	std::cout << pt << std::endl;      
        target_locs.push_back(PMP::locate(pt, mesh));
      } 
    }
  }

  printf("trying to build submesh: \n"); 
  
  auto start = std::chrono::high_resolution_clock::now();
  Triangle_mesh submesh = build_minimum_submesh(source_loc, target_locs, cutoff_rad, mesh); 
  auto end = std::chrono::high_resolution_clock::now();
  auto submesh_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
  std::cout << "time to traverse logic & make submesh: " << submesh_time.count() << "ms" << std::endl;
 
  std::vector<Point_3> source_vec = {source_pt};
  time_big_sequence_tree(submesh, source_vec);
  time_big_sequence_tree(mesh, source_vec);

  std::ofstream out("submesh.off");
  out << submesh << std::endl;

  return 0; 
  
}

