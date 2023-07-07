#ifndef SHIFT_H
#define SHIFT_H


auto normalize(Vector_3 v);
auto getVertexPositions(Triangle_mesh mesh, Triangle_mesh::face_index fIndex);
auto getVertexToRotate(std::vector<Point_3> vs2, std::vector<Point_3> sEdge);
auto rotateAboutSharedAxis(Point_3 target, std::vector<Point_3> axis, double rotAngle);
auto overEdge(Triangle_mesh mesh, Face_location f1, Face_location f2, Point_3 pos, Vector_3 move);
auto shift(Triangle_mesh mesh, Point_3 pos, Vector_3 move);






#endif
