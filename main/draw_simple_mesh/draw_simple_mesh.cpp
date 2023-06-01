#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <fstream>
typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;
int main(int argc, char* argv[])
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("/home/twebb8/projects/cgal/c_curved_space/c_curved_space/main/draw_simple_mesh/heightmap.off");
  Mesh sm;
  
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  CGAL::draw(sm);
  return EXIT_SUCCESS;
}
