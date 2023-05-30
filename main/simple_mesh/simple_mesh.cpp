#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/number_utils.h>
#include <math.h>
#include <fstream>
// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

//height function acting as utility to define the underlying surface
double height(double x, double y) {
  return (sin(x)+sin(y))/2.0;
}

//we can define our surface implicitly via a constraint -- mirroring sphere example
FT heightmap_function(Point_3 p) {
  const FT x = p.x(), y = p.y(), z = p.z();
  return z - height(x,y);
}

int main() {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
  // defining the surface
  Surface_3 surface(heightmap_function,             // pointer to function
                    Sphere_3(CGAL::ORIGIN, 4.)); // bounding sphere
  // Note that "4." above is the *squared* radius of the bounding sphere, which
  // I believe defines the domain of the surface
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     0.1,  // radius bound
                                                     0.1); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
  Surface_mesh sm;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
  std::ofstream out("heightmap.off");
  out << sm << std::endl;
  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
