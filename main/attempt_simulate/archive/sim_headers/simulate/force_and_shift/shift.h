#ifndef SHIFT_H
#define SHIFT_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <math.h>

auto normalize(Vector_3 v);
auto getVertexPositions(Triangle_mesh mesh, Triangle_mesh::face_index fIndex);
auto getVertexToRotate(std::vector<Point_3> vs2, std::vector<Point_3> sEdge);
auto rotateAboutSharedAxis(Point_3 target, std::vector<Point_3> axis, double rotAngle);
auto overEdge(Triangle_mesh mesh, Face_location f1, Face_location f2, Point_3 pos, Vector_3 move);
auto shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move);

#endif
