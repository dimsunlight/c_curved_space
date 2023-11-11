#ifndef UTILS_H
#define UTILS_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
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
typedef Triangle_mesh::Halfedge_index                                   Halfedge_index;
typedef Triangle_mesh::Vertex_index                                     Vertex_index;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Triangle_mesh, FT>                           Face_location;


std::vector<Point_3> getVertexPositions(Triangle_mesh mesh, Face_index fIndex);
std::vector<Vertex_index> getVertexIndices(Triangle_mesh mesh, Face_index fIndex); 
double vectorMagnitude(Vector_3 v);
double vectorMagnitude(Vector_2 v);
Vector_3 normalizer(Vector_3 v);
Vector_2 normalizer(Vector_2 v); 
int findIndex(Point_3 loc,std::vector<Point_3> vec);
auto getSharedElements(std::vector<Point_3> vs1, std::vector<Point_3> vs2);
auto getUnsharedElements(std::vector<Point_3> vertices, std::vector<Point_3> shared);
double angleBetween(Face_index f1, Face_index f2, Triangle_mesh mesh);
std::vector<Point_3> rotateAboutAxis(std::vector<Point_3> targets, std::vector<Point_3> axis, double rotAngle);
Point_3 rotateAboutAxis(Point_3 target, std::vector<Point_3> axis, double rotAngle);
auto createTemporaryMesh(std::vector<Point_3> vertexTrio);
Point_3 getNewXYZLocation(Point_3 flatLocation, Triangle_mesh originalMesh, Triangle_mesh tempMesh, Face_index f2); 
std::vector<Segment_3> createEdgeSegments(std::vector<Point_3> vs);
Face_index getTargetFace(Point_3 pos, Vector_3 toIntersection, Face_index currentSourceFace, Triangle_mesh mesh);
Vector_3 reduceVector(Vector_3 moveVector, double reduceLengthBy);
Vector_3 projectMoveDown(Point_3 source, Vector_3 targetFaceNormal, Vector_3 moveDirection, double remainingDistance);
double triangle_area(Point_3 v1, Point_3 v2, Point_3 v3);
Point_3 project_to_face(std::vector<Point_3> vertices, Point_3 query);

#endif

