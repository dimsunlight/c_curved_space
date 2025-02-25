#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
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

FT torus_function(Point_3 p) {
  //eqn: (c-sqrt(x^2+y^2))^2 + z^2 = a^2
  //for now, defining c and a here so they'll be well-behaved for tests
  double c = 3.0;
  double a = 1.0;
  const FT x2 = p.x()*p.x(), y2 = p.y()*p.y(), z2 = p.z()*p.z();
  return (c-sqrt(x2+y2))*(c-sqrt(x2+y2))+z2 - a*a;
}

int main() {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
  // defining the surface
  Surface_3 surface(torus_function,             // pointer to function
                    Sphere_3(CGAL::ORIGIN, 20.)); // bounding sphere
  // Note that "20." above is the *squared* radius of the bounding sphere!
  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     0.05,  // radius bound
                                                     0.1); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
  Surface_mesh sm;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
  std::ofstream out("torusrb05.off");
  out << sm << std::endl;
  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
