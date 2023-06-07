#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
typedef Kernel::Point_3                                                 Point_3;
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
int main(int argc, char** argv)
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("sims_project/torus.off");
  Triangle_mesh tmesh;
  if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
     !CGAL::is_triangle_mesh(tmesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  Surface_mesh_shortest_path shortest_paths(tmesh);
  // The source point is a 3D point. We must convert it to a face location (that is,
  // a face ID and three barycentric coordinates)
  const Point_3 source_pt(3.51033, 1.91770, 0.000000);
  Face_location source_loc = shortest_paths.locate<AABB_face_graph_traits>(source_pt); // this builds an AABB tree of the mesh
  // Add the source point
  shortest_paths.add_source_point(source_loc.first, source_loc.second);
  // The target will be defined by the first intersection of the following ray with the input mesh
  const Point_3 target_pt(2.08546, 1.66027, 0.942445);
  // If you intend to call `locate` many times, you should cache the AABB tree
  // that is used internally to avoid it being recomputed at each call, as follows:
  AABB_tree tree;
  shortest_paths.build_aabb_tree(tree);
  Face_location target_loc = shortest_paths.locate<AABB_face_graph_traits>(target_pt);

  // Compute the shortest path between the source and the target
  std::vector<Point_3> points;
  shortest_paths.shortest_path_points_to_source_points(target_loc.first, target_loc.second,
                        		  std::back_inserter(points));
  //std::back_inserter above fills points by inserting point by point at the end for each
  
  // Compute the shortest distance between the source and the target
  std::cout << "distance: " << std::endl;
  std::cout << std::get<0>( shortest_paths.shortest_distance_to_source_points(target_loc.first, target_loc.second)) << std::endl;


  // Print the points
  std::cout << points.size() << " ";
  for (std::size_t i = 0; i < points.size(); ++i)
    std::cout << " " << points[i];
  std::cout << std::endl;
  return EXIT_SUCCESS;
  // Print the distance
  
}

/* How the other program creates a polylines file:
   * std::ofstream output("shortest_paths_with_id.polylines.txt");
  vertex_iterator vit, vit_end;
  for ( boost::tie(vit, vit_end) = vertices(tmesh);
        vit != vit_end; ++vit)
  {
    std::vector<Traits::Point_3> points;
    shortest_paths.shortest_path_points_to_source_points(*vit, std::back_inserter(points));
    // print the points
    output << points.size() << " ";
    for (std::size_t i = 0; i < points.size(); ++i)
      output << " " << points[i];
    output << std::endl;
  }*/

