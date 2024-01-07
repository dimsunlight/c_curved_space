#ifndef GET_FORCE_H
#define GET_FORCE_H
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

int time_big_sequence_tree ( const Triangle_mesh &mesh, const std::vector<Point_3> &source);
std::pair<std::vector<double>,std::vector<Vector_3>> calcTangentsAndDistances (
                const Triangle_mesh &mesh, const Point_3 &source, const std::vector<Point_3> &targets, const std::size_t &num_targets);
Vector_3 forceFunction (const float &dist, const Vector_3 &tangent, const double &epsilon, const double &sigma);
Vector_3 force_on_source (const Triangle_mesh &mesh, const Point_3 &source, const std::vector<Point_3> &targets, const std::size_t &num_targets);
Triangle_mesh build_minimum_submesh(const Face_location& source, const std::vector<Face_location>& targets,
                                    const double& cutoff_dist, const Triangle_mesh& mesh);
Vector_3 bounded_region_force(const Triangle_mesh &mesh, const Point_3 &source, const std::vector<Point_3> &targets, 
		              const std::size_t &num_targets, double cutoff_rad);

#endif

