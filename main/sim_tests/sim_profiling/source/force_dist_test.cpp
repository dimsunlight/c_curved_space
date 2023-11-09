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

int main(int argc, char* argv[]) {
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("sims_project/torusrb20.off");
  Triangle_mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  //for (int i = 0; i < 100; i++) std::cout << torus_sample(1,3) << std::endl; //check sampling for accuracy
  
  //now: get two random points on the torus, calculate distance between them, see how long it took to 
  //calculate, and place time to calculate in a bin based on how far apart they were. I think the binning
  //will happen in mathematica -- file output should just be a pair on each line of (dist, time) objects 
  Point_3 point1;
  std::vector<Point_3> point2; 
  std::chrono::duration<long, std::milli> f_time;
  std::pair<std::vector<double>,std::vector<Vector_3>> dist_and_tang; 
  std::size_t targetcount = 1;
  
  const std::string dist_and_time_filename = "dist_and_time.csv";
  std::ofstream dist_and_time_file(dist_and_time_filename);
 
  for (int i = 0; i < 1000; i++) {
    point1 = torus_sample(1,3);
    point2 = {torus_sample(1,3)};

    auto start = std::chrono::high_resolution_clock::now();
    Vector_3 f_on_p = force_on_source(mesh, point1, point2, targetcount);
    auto end = std::chrono::high_resolution_clock::now();
    f_time =  std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    dist_and_tang = calcTangentsAndDistances (mesh, point1, point2, targetcount);
    if (dist_and_time_file.is_open()){
      dist_and_time_file << dist_and_tang.first[0] << ", " <<  f_time.count();
      dist_and_time_file << "\n";
    }
  }

}









