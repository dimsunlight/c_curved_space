#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <iostream>
#include <string>
// default triangulation for Surface_mesher
typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef K::Point_3                                                Point;
typedef K::Vector_3                                               Vector;
typedef CGAL::Surface_mesh<Point>                                 Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor        face_descriptor;

namespace PMP= CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[]) {
  
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");
  Surface_mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }  
  auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>("v:normals",CGAL::NULL_VECTOR).first;
  auto fnormals = mesh.add_property_map<face_descriptor, Vector>("f:normals",CGAL::NULL_VECTOR).first;
  
  PMP::compute_normals (mesh,vnormals,fnormals);

  std::string outputFilename = filename.substr(0, filename.size()-4)+"normals.txt";
  std::ofstream out(outputFilename);
  for(face_descriptor fd: faces(mesh))
    out << fnormals[fd] << "\n";

  out.close();
  return 0;
}
