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

//utility functions
auto normalize(Vector_3 v)
{
  auto const slen = v.x()*v.x() + v.y()*v.y()+v.z()*v.z();
  auto const d = CGAL::approximate_sqrt(slen);
  return v / d;
}

//function definitions
std::pair<std::vector<double>,std::vector<Vector_3>> calcTangentsAndDistances (
		Triangle_mesh mesh, Point_3 source, std::vector<Point_3> targets, std::size_t num_targets) {
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
    Face_location target_loc = shortest_paths.locate<AABB_face_graph_traits>(targets[i],tree);
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
  
  Face_location target_loc = shortest_paths.locate<AABB_face_graph_traits>(target,tree);
  std::vector<Point_3> points;

  shortest_paths.shortest_path_points_to_source_points(target_loc.first, target_loc.second, std::back_inserter(points));
  
  double distance;
  Vector_3 tangent;
  
  distance = std::get<0>(shortest_paths.shortest_distance_to_source_points(target_loc.first,target_loc.second));
  tangent = normalize(Vector_3(points[points.size()-2],points[points.size()-1]));

  return std::make_pair(distance,tangent);
}

auto forceFunction (float dist, Vector_3 tangent, double epsilon, double sigma) { 
  //once again, a=1 and c=3

  Vector_3 normalizedTangent = normalize(tangent);

  //lennard-jones potential is 
  //4 epsilon((sigma/r)^12 - (sigma/r)^6) 
  //associated force is 
  //4 epsilon (12 sigma^12 / r^13 - 6 sigma ^6 / r^7)nabla(r)
  //where r is the distance function on the surface 
  
  Vector_3 force = 4*epsilon*(12*pow(sigma,12)/pow(dist,13)-6*pow(sigma,6)/pow(dist,7))*normalizedTangent;

  return force;
}


Vector_3 force_on_source (Triangle_mesh mesh, Point_3 source, std::vector<Point_3> targets, std::size_t num_targets) {
  //create list of distances and path tangents between the source particle and the targets
  std::pair<std::vector<double>, std::vector<Vector_3>> distancesAndTangents = calcTangentsAndDistances(mesh, source, targets, num_targets);
  std::vector<double> distances = distancesAndTangents.first;
  std::vector<Vector_3> tangents = distancesAndTangents.second;
  //we could either do this loop within forceFunction or here -- shouldn't be hugely different

  //L-J parameters
  double epsilon = 1;
  double sigma = .8;
  Vector_3 force= Vector_3(0,0,0); //initialize to zero to avoid redefinition --
                                   //also handles case of no neighbors

  for (std::size_t i = 0; i < distances.size(); i++) {
    force+= forceFunction(distances[i],tangents[i], epsilon, sigma);
  }

  return force;
}


std::vector<std::pair<Point_3,Point_3>> pathEndpointDatabase () {
  
  std::vector<std::pair<Point_3,Point_3>> pathEndpoints;	
  pathEndpoints = {std::make_pair(Point_3(4.,0.,0.),Point_3(2.16121,3.36588,0.)),
  std::make_pair(Point_3(3.98769,0.,0.156434),Point_3(2.17237,3.35625,-0.064079)),
  std::make_pair(Point_3(3.95106,0.,0.309017),Point_3(2.20543,3.32834,-0.120505)),
  std::make_pair(Point_3(3.89101,0.,0.453990),Point_3(2.25906,3.28499,-0.162016)),
  std::make_pair(Point_3(3.80902,0.,0.587785),Point_3(2.33067,3.23026,-0.182028)),
  std::make_pair(Point_3(3.70711,0.,0.707107),Point_3(2.41598,3.16861,-0.174817)),
  std::make_pair(Point_3(3.58779,0.,0.809017),Point_3(2.50848,3.10381,-0.135697)),
  std::make_pair(Point_3(3.45399,0.,0.891007),Point_3(2.59900,3.03811,-0.0614595)),
  std::make_pair(Point_3(3.30902,0.,0.951057),Point_3(2.67585,2.97159,0.0487022)),
  std::make_pair(Point_3(3.15643,0.,0.987688),Point_3(2.72561,2.9022,0.191874)),
  std::make_pair(Point_3(3.,0.,1.)           ,Point_3(2.73487,2.82632,0.360174)),
  std::make_pair(Point_3(2.84357,0.,0.987688),Point_3(2.69261,2.73984,0.540303)),
  std::make_pair(Point_3(2.69098,0.,0.951057),Point_3(2.59273,2.63931,0.714381)),
  std::make_pair(Point_3(2.54601,0.,0.891007),Point_3(2.43577,2.52291,0.862032)),
  std::make_pair(Point_3(2.41221,0.,0.809017),Point_3(2.2295,2.39104,0.963081)),
  std::make_pair(Point_3(2.29289,0.,0.707107),Point_3(1.98836,2.24654,1.)),
  std::make_pair(Point_3(2.19098,0.,0.587785),Point_3(1.73217,2.09499,0.959514)),
  std::make_pair(Point_3(2.10899,0.,0.45399),Point_3(1.48513,1.9456,0.833607)),
  std::make_pair(Point_3(2.04894,0.,0.309017),Point_3(1.27512,1.81296,0.621361)),
  std::make_pair(Point_3(2.01231,0.,0.156434),Point_3(1.13184,1.71801,0.333723)),
  std::make_pair(Point_3(2.,0.,0.),Point_3(1.0806,1.68294,0.))};

  return pathEndpoints;
}

std::vector<Vector_3> pathEndpointTangents() {
  std::vector<Vector_3> tangents;
  tangents = {Vector_3(0., 2., 0.), 
	      Vector_3(0., 1.99384, 0.),
	      Vector_3(0., 1.97553, 0.), 
	      Vector_3(0., 1.9455, 0.),
	      Vector_3(0., 1.90451, 0.),
	      Vector_3(0., 1.85355, 0.),
	      Vector_3(0., 1.79389, 0.),
	      Vector_3(0., 1.72699, 0.),
	      Vector_3(0., 1.65451, 0.),
	      Vector_3(0., 1.57822, 0.),
	      Vector_3(0., 1.5, 0.),
	      Vector_3(0., 1.42178, 0.),
	      Vector_3(0., 1.34549, 0.),
	      Vector_3(0., 1.27301, 0.),
	      Vector_3(0., 1.20611, 0.),
	      Vector_3(0., 1.14645, 0.),
	      Vector_3(0., 1.09549, 0.),
	      Vector_3(0., 1.0545,  0.),
	      Vector_3(0., 1.02447, 0.),
	      Vector_3(0., 1.00616, 0.),
	      Vector_3(0., 1., 0.)};
  
  return tangents;
}

std::vector<double> pathDistances() { 
  std::vector<double> distances= {4.0,3.98769,3.95106,3.891,3.80901,
	       3.70711,3.58778,3.45399,3.30902,
	       3.15644,3.00000,2.84357,2.69098,
	       2.54601,2.41222,2.29290,2.19099,
	       2.10900,2.04894,2.01231,2.0000};

  return distances;
}

int main(int argc, char* argv[]) {
 //prepare path endpoints and source tangents
 std::vector<std::pair<Point_3,Point_3>> endpoint_db = pathEndpointDatabase();
 std::vector<Vector_3> tangents = pathEndpointTangents();
 std::vector<double> distances = pathDistances(); 
 
 //get input mesh from command line argument
 
 const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("sims_project/heightmap.off");
 Triangle_mesh tmesh;
 if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
   !CGAL::is_triangle_mesh(tmesh))
 {
   std::cerr << "Invalid input file." << std::endl;
   return EXIT_FAILURE;
 }
 std::cout << "mesh loaded, number of vertices " << tmesh.number_of_vertices() << std::endl; 

 //test time to calculation for every path in the database
 double trueDistance;
 double calcedDistance;
 double tangentOverlap;
 std::vector<double> path_calc_times;
 std::vector<std::pair<double,double>> distance_and_tangent_error;
 clock_t t;
 std::pair<double,Vector_3> current_d_and_t;
 std::size_t counter = 0;
 for (std::pair<Point_3,Point_3> st_pair: endpoint_db) {
   t = clock();
   current_d_and_t = calcTangentandDistance(tmesh,st_pair.first,st_pair.second);
   clock_t end = clock();
   path_calc_times.push_back( ((double) (clock() - t)) / CLOCKS_PER_SEC);
   calcedDistance = current_d_and_t.first;
   trueDistance = distances[counter];
   tangentOverlap = CGAL::scalar_product(current_d_and_t.second,normalize(tangents[counter]));
   distance_and_tangent_error.push_back(std::make_pair(trueDistance-calcedDistance, tangentOverlap));
   counter++;
 }
 
 double averageTangentOverlap = 0;
 double averageDistanceError = 0;
 double averageTimeToCalc = 0;
 counter = 0;//reusing old counter so the variable name is clear
 //std::cout << "LISTING INDIVIDUAL ERRORS & TIMES:" << std::endl; 
 for (std::pair<double,double> errors: distance_and_tangent_error) {
   //std::cout << "distance error " << errors.first << ", tangent overlap " << abs(errors.second) << std::endl;
   //std::cout << "calculation time " << path_calc_times[counter] << std::endl;
   averageDistanceError += errors.first;
   averageTangentOverlap += abs(errors.second);
   averageTimeToCalc += path_calc_times[counter]; 
   counter++;  
 }
 std::cout << "-------------;" << std::endl;
 averageDistanceError = averageDistanceError/(counter);
 averageTangentOverlap = averageTangentOverlap/(counter);
 averageTimeToCalc = averageTimeToCalc/(counter);

 std::cout << "average distance error " << averageDistanceError << std::endl;
 std::cout << "average tangent overlap " << averageTangentOverlap << std::endl;
 std::cout << "average calculation time: " << averageTimeToCalc << std::endl;
 

 /*
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
