#ifndef SHIFT_H
#define SHIFT_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
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
#include <math.h> //including this to do simple cos and sine, but maybe imprecise compared to cgal-originating routines
//types and names
typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::FT                                                           FT;
typedef K::Point_2                                                      Point_2;
typedef K::Ray_2                                                        Ray_2;
typedef K::Point_3                                                      Point_3;
typedef K::Vector_3                                                     Vector_3;
typedef K::Ray_3                                                        Ray_3;
typedef K::Segment_3                                                    Segment_3;
typedef K::Intersect_3                                                  Intersect_3;
typedef K::Triangle_3                                                   Triangle_3;
typedef K::Segment_2                                                    Segment_2;
typedef K::Point_2                                                      Point_2;
typedef K::Vector_2                                                     Vector_2;
typedef CGAL::Surface_mesh<Point_3>                                     Triangle_mesh;
typedef Triangle_mesh::Face_index                                       Face_index;
typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor  vertex_descriptor;
typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor    face_descriptor;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;

auto normalize(Vector_3 v);
auto getVertexPositions(Triangle_mesh mesh, Triangle_mesh::face_index fIndex);
auto getVertexToRotate(std::vector<Point_3> vs2, std::vector<Point_3> sEdge);
auto rotateAboutAxis(Point_3 target, std::vector<Point_3> axis, double rotAngle);
Point_3 getNewXYZLocation(Point_3 flatLocation, Triangle_mesh originalMesh, Triangle_mesh tempMesh, face_descriptor f2);
Vector_3 projectMoveDown(Point_3 source, Vector_3 targetFaceNormal, Vector_3 moveDirection, double remainingDistance); 
double triangle_area(Point_3 v1, Point_3 v2, Point_3 v3);
Point_3 project_to_face(std::vector<Point_3> vertices, Point_3 query);
Point_3 find_intersection_baryroutine(Point_3 source, Point_3 target,  std::vector<Point_3> faceVertices); 
Point_3 shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move);


#endif

