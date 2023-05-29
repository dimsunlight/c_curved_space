#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <cmath>
typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

double height(double x, double y) {
  return (sin(x)+sin(y))/2.0;
}

double[10] listOfPoints() {
  return 1;
}



int main()
{
  Mesh m;
  // Add the points as vertices
  vertex_descriptor u = m.add_vertex(K::Point_3(0,1,0));
  vertex_descriptor v = m.add_vertex(K::Point_3(0,0,0));
  vertex_descriptor w = m.add_vertex(K::Point_3(1,1,0));
  vertex_descriptor x = m.add_vertex(K::Point_3(1,0,0));
  m.add_face(u,v,w);
  face_descriptor f = m.add_face(u,v,x);
  if(f == Mesh::null_face())
  {
    std::cerr<<"The face could not be added because of an orientation error."<<std::endl;
    f = m.add_face(u,x,v);
    assert(f != Mesh::null_face());
  }
  return 0;
}
