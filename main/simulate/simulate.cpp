/* Author: Toler H. Webb
 * Description: CGAL-based routine for simulating the interaction of a 
 * single particle pair. 
 */

//includes
#include <CGAL>
#include <get_force.h>
#include "shift.h"
#include "get_force.h"

int main () {
  //have defaults for both loading functionalities below
  Triangle_mesh mesh = loadMesh(input_mesh);
  Point_3 particle_locations[];
  particle_locations = create_particles_from_xyz(input_xyz); 
  Point_3 location_buffer[particle_locations.size()];

  std::vector<std::pair<Point_3,std::vector<Point_3>>> particles_with_neighbors = get_neighbors;
  
  //can make this multi-step simply by enclosing into timestep loop
  std::pair<Point_3,std::vector<Point_3>> particle_and_neighbors;
  Vector_3 f_on_p;
  for (std::size_t i = 0; i < location_buffer.size();i++) {
    particle_and_neighbors = particles_with_neighbors[i];
    f_on_p = force_on_source(mesh,particle_and_neighbors.first,particle_and_neighbors.second,
		    particle_and_neighbors.second.size());
    locationBuffer[i] = shift(mesh, particle_and_neighbors.first,f_on_p);
  }
  //here, we could store location_buffers in a std::vector to use them
  //for animation later
  particle_locations = location_buffer;
  

  return 0;
}





