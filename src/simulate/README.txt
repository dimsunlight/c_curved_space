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
