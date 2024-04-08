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
typedef std::mt19937                                                    MTgenerator; //mersenne twister generator
typedef std::uniform_real_distribution<float>                           distribution;

int main(int argc, char* argv[]) {
    const std::string read_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("torus_isotropic.off");
    Triangle_mesh mesh;
    if(!CGAL::IO::read_polygon_mesh(read_filename,mesh) || !CGAL::is_triangle_mesh(mesh))
    {
      std::cerr << "Invalid input file." << std::endl;
      return EXIT_FAILURE;
    }
    Point_3 p1 = Point_3(3.52,1.15,0.71);
    Point_3 p2 = Point_3(1.08, 3.33, 0.87);
   
    Surface_mesh_shortest_path shortest_paths(mesh);
    AABB_tree tree;
    shortest_paths.build_aabb_tree(tree);
    
    Face_location loc1 = shortest_paths.locate<AABB_face_graph_traits>(p1,tree);
    Face_location loc2 = shortest_paths.locate<AABB_face_graph_traits>(p2,tree);

    std::cout << "face locations" << std::endl;
    std::cout << loc1.first << ", " << loc1.second[0] << " " << loc1.second[1] << " " << loc1.second[2] << std::endl;
    std::cout << loc2.first << ", " << loc2.second[0] << " " << loc2.second[1] << " " << loc2.second[2] << std::endl;

    shortest_paths.add_source_point(loc1);

    shortest_paths.build_sequence_tree();

    std::vector<Point_3> points;
    shortest_paths.shortest_path_points_to_source_points(loc2.first, loc2.second, std::back_inserter(points)); 

    std::ofstream ps("example_geodesic.txt");
    ps << "{";  
    for (Point_3 p: points) { 
        ps << "{" << p.x() << ", " << p.y() << ", " << p.z() << "}, "; 
    }
    ps << "}"; 
    return 0; 

}
