# Curved Space Molecular Dynamics simulations using C++ 
(Repo is very much a work in progress -- there's lots of junk to be removed and helper code to be written before it's ready for public consumption) 

This repo stores code meant to execute curved space simulations primarily with the help of CGAL. 

##**Main**

Main is a replica of the old filestructure, in case I've made any mistakes reorganizing. It has both copies of the
code I'm currently working on and deprecated files from when I was first writing the repository. I'll know where to find things in there, but there's no reason anyone looking at the repository should need to look at it. 

##**ex** 

Placeholder folder for examples, once I write some. 

##**src**

Main executable folder for the most recent stable version of the code. 
  - classes is a (mostly) empty folder for when/if I restructure the code to be object-oriented (and thereby more
  modular so I can include more built-in functionality). 
  - mesher has several functions built in for mesh generation & related functionality. I'll be cleaning up the version
  in src later to remove extraneous files -- e.g. sphere_triangulation_example, simple_mesh_plan, or compute_normals
  should proabbly exist in test. 
  - meshes is a collection of meshes I've generated to make it easy to stay consistent across simulation runs for tests.
  - simulate is the current stable simulation code, with the most up-to-date shift and force functions. As I write
  this, it has testing files in it. I'll clean up the src version when I clean up the mesher so extraneous files are 
  deleted. 

##**test** 

collection of test files. Nothing here is guaranteed to be stable. 
- cgal_distance is a record of my toying with CGAL's Xin and Wang distance algorithm. 
- force has a functional isolated version of the force function, with testing-focused code in its main. 
- shift is the equivalent isolated, up-to-date shift function (as of now). That shift function executes a path-bending routine for finding geodesics on the mesh. The folder within shift, old, has some previous attempts to write shift correctly. I've found them useful when debugging. 
- first_sim has a lot of deprecated code in it from when I was writing my first simulation functions. 
- simulate has the current, testing/timing focused versions of the simulation code. At present this is not 
different from what's in src. 
- timing_pieces has its own ReadMe. It's meant to provide a bunch of functions to time and test the subordinate 
components of the simulation routine. 

##**How you should use this repo** 

If you're curious about my work (or if you're my advisor), a CGAL installation should be enough to let you run
the simulation routine in src. I'll put some instructions there & in test, because it's straightforward to repeat but nontrivial. Otherwise, I'll indicate here parts of test that might need attention. 


## TODO: 
- Clean up src/simulate's files so that it's only the stable simulate that can run a simulation on a mesh. 
- Clean up simulate files within test. 
- Fix MT sampling within random_pos_simulate (within test file -- unnecessary, for now, in src) so that we get 
sample that's actually random. 
- Create some examples. 
- Make code more user friendly (not having to vim into force function to change what force we use, for example, or
vim into simulate function to change timesteps)
