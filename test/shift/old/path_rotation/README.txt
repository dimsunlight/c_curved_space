This code represents the first set of path rotation shift function implementations. 
The goal of this was to rotate paths into the faces that they traveled into over the course
of a shift. In principle, this should produce a geodesic. 

A couple flaws here: 
-get_target_faces is wrong. It should use topological info of the mesh to decide what face to travel into. 
-utility function structure is confusing & disorganized
-intersection tests are of varying quality -- best to just use barycentric coordinates to make it easy