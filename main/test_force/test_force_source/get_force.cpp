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

//function definitions
std::pair<std::vector<double>,std::vector<Vector_3>> calcTangentsAndDistances (
		Triangle_mesh mesh, Point_3 source, Point_3 targets[], std::size_t num_targets) {
  Surface_mesh_shortest_path shortest_paths(mesh);
  AABB_tree tree;
  shortest_paths.build_aabb_tree(tree);
  //convert source point to barycentric coordinates via locate
  const Point_3 source_pt = source;
  Face_location source_loc = shortest_paths.locate<AABB_face_graph_traits>(source_pt);
  shortest_paths.add_source_point(source_loc.first,source_loc.second);

  std::vector<double> distances;
  std::vector<Vector_3> tangents; 
  for (std::size_t i = 0; i < num_targets; i++) {
    Face_location target_loc = shortest_paths.locate<AABB_face_graph_traits>(targets[i]);
    std::vector<Point_3> points; 
    shortest_paths.shortest_path_points_to_source_points(target_loc.first, target_loc.second, std::back_inserter(points));
    distances.push_back(std::get<0>( shortest_paths.shortest_distance_to_source_points(target_loc.first, target_loc.second)));
    //path goes from target to source -- so if we want to know path tangent at 
    //source for force calculation, we must use the *end* of points[]
    tangents.push_back(Vector_3(points[points.size()-2],points[points.size()-1]));
  }
 
  return std::make_pair(distances,tangents);
}

std::pair<double,Vector_3> calcTangentandDistance (Triangle_mesh mesh, Point_3 source, Point_3 target) {
  
  Surface_mesh_shortest_path shortest_paths(mesh);

  const Point_3 source_pt = source;
  Face_location source_loc = shortest_paths.locate<AABB_face_graph_traits>(source_pt);
  shortest_paths.add_source_point(source_loc.first,source_loc.second);
  
  AABB_tree tree;
  shortest_paths.build_aabb_tree(tree);
  
  Face_location target_loc = shortest_paths.locate<AABB_face_graph_traits>(target);
  std::vector<Point_3> points;

  shortest_paths.shortest_path_points_to_source_points(target_loc.first, target_loc.second, std::back_inserter(points));
  
  double distance;
  Vector_3 tangent;
  
  std::cout << "interpolation point 1: " << points[points.size()-2] << std::endl;
  std::cout << "interpolation point 2: " << points[points.size()-1] << std::endl;
  distance = std::get<0>(shortest_paths.shortest_distance_to_source_points(target_loc.first,target_loc.second));
  tangent = Vector_3(points[points.size()-2],points[points.size()-1]);

  return std::make_pair(distance,tangent);
}

auto normalize(Vector_3 v)
{
  auto const slen = v.x()*v.x() + v.y()*v.y()+v.z()*v.z();
  auto const d = CGAL::approximate_sqrt(slen);
  return v / d;
}

auto forceFunction (float dist, Vector_3 tangent) { 
  //once again, a=1 and c=3
  const float epsilon = 1;
  const float sigma = 1; //i know, large particles
  Vector_3 normalizedTangent = normalize(tangent);

  //lennard-jones potential is 
  //4 epsilon((sigma/r)^12 - (sigma/r)^6) 
  //associated force is 
  //4 epsilon (12 sigma^12 / r^13 - 6 sigma ^6 / r^7)nabla(r)
  //where r is the distance function on the surface 
  
  Vector_3 force = 4*epsilon*(12*pow(sigma,12)/pow(dist,13)-6*pow(sigma,6)/pow(dist,7))*normalizedTangent;

  return force;
}


struct pointdatabase {
  Point_3 array[21][2];
};

struct pointdatabase pathEndpointDatabase () {
  
  struct pointdatabase pathEndpoints;	
  pathEndpoints.array = {{Point_3(4.,0.,0.),Point_3(2.16121,3.36588,0.)},
  {Point_3(3.98769,0.,0.156434),Point_3(2.17237,3.35625,-0.064079)},
  {Point_3(3.95106,0.,0.309017),Point_3(2.20543,3.32834,-0.120505)},
  {Point_3(3.89101,0.,0.453990),Point_3(2.25906,3.28499,-0.162016)},
  {Point_3(3.80902,0.,0.587785),Point_3(2.33067,3.23026,-0.182028)},
  {Point_3(3.70711,0.,0.707107),Point_3(2.41598,3.16861,-0.174817)},
  {Point_3(3.58779,0.,0.809017),Point_3(2.50848,3.10381,-0.135697)},
  {Point_3(3.45399,0.,0.891007),Point_3(2.59900,3.03811,-0.0614595)},
  {Point_3(3.30902,0.,0.951057),Point_3(2.67585,2.97159,0.0487022)},
  {Point_3(3.15643,0.,0.987688),Point_3(2.72561,2.9022,0.191874)},
  {Point_3(3.,0.,1.)           ,Point_3(2.73487,2.82632,0.360174)},
  {Point_3(2.84357,0.,0.987688),Point_3(2.69261,2.73984,0.540303)},
  {Point_3(2.69098,0.,0.951057),Point_3(2.59273,2.63931,0.714381)},
  {Point_3(2.54601,0.,0.891007),Point_3(2.43577,2.52291,0.862032)},
  {Point_3(2.41221,0.,0.809017),Point_3(2.2295,2.39104,0.963081)},
  {Point_3(2.29289,0.,0.707107),Point_3(1.98836,2.24654,1.)},
  {Point_3(2.19098,0.,0.587785),Point_3(1.73217,2.09499,0.959514)},
  {Point_3(2.10899,0.,0.45399),Point_3(1.48513,1.9456,0.833607)},
  {Point_3(2.04894,0.,0.309017),Point_3(1.27512,1.81296,0.621361)},
  {Point_3(2.01231,0.,0.156434),Point_3(1.13184,1.71801,0.333723)},
  {Point_3(2.,0.,0.),Point_3(1.0806,1.68294,0.)}};



  return pathEndpoints;
}


int main(int argc, char* argv[]) {
 struct pointdatabase endpoints;
 endpoints = pathEndpointDatabase();

 std::cout << endpoints[1][1] << std::endl;
 std::cout << endpoints[7][2] << std::endl;
 std::cout << endpoints[20][2] <<std::endl;
 //get input mesh from command line argument
 /*
 const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("sims_project/heightmap.off");
 Triangle_mesh tmesh;
 if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
   !CGAL::is_triangle_mesh(tmesh))
 {
   std::cerr << "Invalid input file." << std::endl;
   return EXIT_FAILURE;
 }

 clock_t t;

 Point_3 source_pt(3.51033,1.91770,0.000000);
 Point_3 target_pt(2.08546,1.66027,0.942445);
 Point_3 target_pts[] = {target_pt};

 std::cout << "with single-pair function" << std::endl;

 t=clock();
 std::pair<double, Vector_3> oneCalc = calcTangentandDistance (tmesh, source_pt, target_pt);
 clock_t timeTaken = ( clock() -t);

 std::cout << "Distance:" << oneCalc.first << std::endl;
 std::cout << "Force Grad: " << oneCalc.second << std::endl;
 std::cout << "Clocks per second: " << CLOCKS_PER_SEC << std::endl;;
 std::cout << "Time taken: " << timeTaken << " clicks and "  << timeTaken/CLOCKS_PER_SEC << " seconds." << std::endl;

 Vector_3 force = forceFunction(oneCalc.first, oneCalc.second);
 std::cout << "The force is " << force << std::endl;
 */
 return 0;
}
