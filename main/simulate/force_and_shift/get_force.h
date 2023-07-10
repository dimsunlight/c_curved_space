#ifndef GET_FORCE_H
#define GET_FORCE_H

std::pair<std::vector<double>,std::vector<Vector_3>> calcTangentsAndDistances (
                Triangle_mesh mesh, Point_3 source, Point_3 targets[], std::size_t num_targets);
auto forceFunction (float dist, Vector_3 tangent, double epsilon, double sigma);
auto force_on_source (Triangle_mesh mesh, Point_3 source, Point_3 targets[], std::size_t num_targets);

#endif
