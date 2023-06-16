//Author: Toler H. Webb
//Description: CGAL code to calculate the force between two particles at given locations. 
//             Particle shape information is stored in the force function. 

//includes:
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

//typedefs:
typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
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
calcTangentsAndDistances (Surface_mesh mesh, Point_3 source, Point_3 targets[]) {
  Surface_mesh_shortest_path shortest_paths(mesh);
  AABB_tree tree;
  shortest_paths.build_aabb_tree(tree);
  std::vector<double> distances;
  std::vector<Vector_3>[targets.size()] tangents;
  for (std::size_t i = 0; i < targets.size(); ++i) {
    Face_location target_loc = shortest_paths.locate<AABB_face_graph_traits>(targets[i]);
    std::vector<Point_3> points; 
    shortest_paths.shortest_path_points_to_source_points(target_loc.first, 
		    target_loc.second, std::back_inserter(points));
    distances.push_back(std::get<0>( shortest_paths.shortest_distance_to_source_points(
			    target_loc.first, target_loc.second));
    //path goes from target to source -- so if we want to know path tangent at 
    //source for force calculationi, we must use the *end* of points[]
    tangents.push_back(Vector_3(points[points.size()-2],points[points.size()-1])
  }
 
  return std::make_pair(distances,tangents);
}

calcTangents (path) {
  return pair(tangent1,tangent2);
}

forceFunction (int argc, char* argv[]) {
  return force;
}


int main(int argc, char* argv[]) {
 //get input mesh from command line argument
 Surface_mesh mesh;
 if(!PMP::IO::read_polygon_mesh(filename, mesh))
 {
   std::cerr << "Invalid input." << std::endl;
   return 1;
 }

 Point_3 source_pt(3.51033,1.91770,0.000000);
 Point_3 target_pt(2.08546,1.66027,0.942445);

 calcPath(point1,point2)
	
 return 0;
}
