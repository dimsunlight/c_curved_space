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


std::pair<std::vector<double>,std::vector<Vector_3>> calcTangentsAndDistances (
                Triangle_mesh mesh, Point_3 source, Point_3 targets[], std::size_t num_targets);
auto forceFunction (float dist, Vector_3 tangent, double epsilon, double sigma);
auto force_on_source (Triangle_mesh mesh, Point_3 source, Point_3 targets[], std::size_t num_targets);

#endif
