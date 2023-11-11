This suite of files is designed to evaluate the accuracy of my force & shift
generating functions and their execution time, individually.


To test the effects of curvature on force, I've created a database of 
paths of varying curvature on the torus. The following list of 
geodesic paths is defined by endpoint; each has length 2 on the smooth
torus, so we can always evaluate the error in the distance calculation
by comparing to that true value. 
{
{{4.,0.,0.},{2.16121,3.36588,0.}},
{{3.98769,0.,0.156434},{2.17237,3.35625,-0.064079}},
{{3.95106,0.,0.309017},{2.20543,3.32834,-0.120505}},
{{3.89101,0.,0.45399},{2.25906,3.28499,-0.162016}},
{{3.80902,0.,0.587785},{2.33067,3.23026,-0.182028}},
{{3.70711,0.,0.707107},{2.41598,3.16861,-0.174817}},
{{3.58779,0.,0.809017},{2.50848,3.10381,-0.135697}},
{{3.45399,0.,0.891007},{2.599,3.03811,-0.0614595}},
{{3.30902,0.,0.951057},{2.67585,2.97159,0.0487022}},
{{3.15643,0.,0.987688},{2.72561,2.9022,0.191874}},
{{3.,0.,1.},{2.73487,2.82632,0.360174}},
{{2.84357,0.,0.987688},{2.69261,2.73984,0.540303}},
{{2.69098,0.,0.951057},{2.59273,2.63931,0.714381}},
{{2.54601,0.,0.891007},{2.43577,2.52291,0.862032}},
{{2.41221,0.,0.809017},{2.2295,2.39104,0.963081}},
{{2.29289,0.,0.707107},{1.98836,2.24654,1.}},
{{2.19098,0.,0.587785},{1.73217,2.09499,0.959514}},
{{2.10899,0.,0.45399},{1.48513,1.9456,0.833607}},
{{2.04894,0.,0.309017},{1.27512,1.81296,0.621361}},
{{2.01231,0.,0.156434},{1.13184,1.71801,0.333723}},
{{2.,0.,1.22465*10^-16},{1.0806,1.68294,1.22465*10^-16}}
}

The database of meshes, included in the build directory, 
is organized by ``accuracy score''. The score is a measure of
how well the local face normals of each mesh correspond 
to the true normal of the smooth torus in the vicinity 
of that face. 

I test both the time-to-completion of distance 
calculation on meshes with various levels of refinement
and the accuracy of that distance calculation; ultimately, 
I aim to have a good sense of how to control mesh quality
and manage the sacrifices involved with using coarser/finer
meshes. 

