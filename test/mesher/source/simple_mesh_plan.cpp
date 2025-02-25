#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3; //three-dimensional point structure
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

//height function acting as utility to define the underlying surface
double height(double x, double y) {
  return (sin(x)+sin(y))/2.0;
}

double[10] Generate() {
  return 1;
}

void createPoints() {
  return 0;
}

void pointsToMesh() {
  return 0; 
}

int main()
{
 double points[10] = GenerateTenPoints();
 Mesh theMesh;
 theMesh = makePointsIntoMesh(points);
 return theMesh;
}
