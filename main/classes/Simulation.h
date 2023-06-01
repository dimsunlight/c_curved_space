//Author: Toler H. Webb
//Description: Class object to store simulation details. Functions stored elsewhere. 

#ifndef Simulation_H
#define Simulation_H
#include<CGAL::CGAL>

using namespace std;

typedef CGAL::Surface_mesh<Point_3> Surface_mesh;


class Simulation{
  public:
    // accessors
    std::vector  xyz_positions() { return xyzPositions_; };
    std::vector  bary_positions() { return baryPositions_; };

    // mutators
    void define_mesh(Surface_mesh inMesh) { underlying_mesh_ = inMesh; };
    void set_positions(double[] cartPos) { xyzPositions_ = cartPos; };

    //class functions
    void generateForces(energy_function);
    void generateDisplacements(forces);
    void updatePositions(dr);
    void convertToBary();
    void convertToXYZ();
    
  private: 
    Surface_mesh underlying_mesh;
    double[] baryPositions;
    double[] xyzPositions;
}
