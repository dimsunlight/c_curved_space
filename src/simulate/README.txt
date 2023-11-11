I will assume that you've got CGAL installed from source (rather than with, e.g.,
sudo apt-get). This code works with CGAL 5.6, specifically. 

To use this simulate: 

-go to build. 
-remove whatever CMakeCache.txt file I've left in there while running my tests. 
-type cmake ../source to link together the build and source folders. 
-type cmake-gui ../source to open the cmake gui. 
-Where it says the CGAL Directory is NOT FOUND, click to go into your file browser. 
-Select the folder where CGAL is installed (so the one called, for example, CGAL 5.6). Make sure that 
folder is specified as the path to CGAL. 
-Click configure until all the red goes away, then click generate. 
-Close out of the gui. 
-type cmake -DCMAKE_BUILD_TYPE="Release" ../source (if you don't configure cmake this way, your simulation will potentially not finish in your lifetime -- it defaults to debug mode) 
-type make
-run executables as you like. For running a simulation, you might use the command line argument:
./simulate "torusrb20.off" "ten_particles.xyz" to run simulate on a coarse torus with the sample particle configuration
in ten_particles.xyz. 

At present, compilation is really slow (and that's only an issue because I have a short attention span). I'm 
very open to tips that would speed up compilation! 
