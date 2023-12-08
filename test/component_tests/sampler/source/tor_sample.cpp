#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <math.h>
#include <chrono>
#include <ratio>
#include <thread>
#include <random>
#include "utils.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     KernelWithSqrt;
typedef Kernel::Point_3                                                 Point_3;
typedef Kernel::Vector_3                                                Vector_3;
typedef std::mt19937                                                    MTgenerator; //mersenne twister generator
typedef std::uniform_real_distribution<float>                           distribution;

MTgenerator make_generator(int seed)
{
  //create random number generator using seed parameter. 

  std::random_device                  rand_dev;
  MTgenerator                         generator(rand_dev());
  generator.seed(seed); //seeded for reproducibility -- make sure we use the same sequence every time :)  

  return generator;
}

distribution distribution_of_range(float range_min, float range_max) {
  distribution distro(range_min, range_max);
  return distro; 
}

float sinplane_weight(float x, float y) {
  return 1/sqrt(1+cos(x)*cos(x)/4.0 + cos(y)*cos(y)/4.0);
}

float torus_weight(float a, float c, float u, float v) {
  //get torus weight for sampling against [0,1] distribution --
  //so we divide by radius to compare to ``unit'' torus
  return (c+a*cos(v))/(c+a);
}

Point_3 sinplane_sample(MTgenerator& gen) {
  //pass by address above to ensure that we sequentially sample 
  //from the generator without using the same starting point.
  float sp_ext = 5.0;
  distribution SPdistro = distribution_of_range(-sp_ext,sp_ext);
  distribution weightDistro = distribution_of_range(0.0,1.0); 
  float X;
  float Y;
  float W;
  float weight;
  Point_3 returnPoint = Point_3(NULL,NULL,NULL);
  while (returnPoint == Point_3(NULL,NULL,NULL)) {
    //goal -- reference gen by address so each call to gen advances the 
    //rng forward once
    X = SPdistro(gen);
    Y = SPdistro(gen);
    W = weightDistro(gen);
    weight = sinplane_weight(X,Y);
    if (W <= weight) {
      returnPoint = Point_3(X,Y,sin(X));
      return returnPoint;
    }
  }
  return returnPoint;
}

Point_3 torus_sample(float a, float c, MTgenerator& gen) {
  //rejection method sample on a torus; that is, reject anything above the corresponding weight. 
  //we pass a mersenne-twister generator by reference 
  distribution angularDistro = distribution_of_range(0.0, 2*M_PI);
  distribution weightDistro = distribution_of_range(0.0,1.0); 
  float U;
  float V;
  float W;
  float weight;
  Point_3 returnPoint = Point_3(NULL,NULL,NULL);
  while (returnPoint == Point_3(NULL,NULL,NULL)) {
    U = angularDistro(gen);
    V = angularDistro(gen);
    W = weightDistro(gen);
    weight = torus_weight(a,c,U,V);
    if (W <= weight) {
      returnPoint = Point_3((c+a*cos(V))*cos(U), (c+a*cos(V))*sin(U), a*sin(V));
      return returnPoint;
    }
  }
  //this should never happen, but just to avoid false error messages
  return Point_3(NULL,NULL,NULL);//filler point to check if we have sampled at all
}

std::vector<Point_3> n_torus_sample_points(std::size_t n, float a, float c, int seed=1997) {
  //I use std::vectors of point_3 objects for positions, so we loop
  //to create a similar vector for n random positions
  std::vector<Point_3> sample_points = {};
  sample_points.reserve(n);
  MTgenerator gen = make_generator(seed); 
  
  for (std::size_t i=0; i < n; i++) {
    sample_points.push_back(torus_sample(a,c, gen));
  }
  return sample_points;  
}

std::vector<Point_3> n_sinplane_sample_points(std::size_t n, int seed = 1997) {
  std::vector<Point_3> sample_points = {};
  sample_points.reserve(n);
  MTgenerator gen = make_generator(seed); 

  for (std::size_t i=0; i < n; i++) {
    sample_points.push_back(sinplane_sample(gen));
  }
  return sample_points;
}

int main(int argc, char* argv[]) {
 
  std::string filename = "tor_500pts.csv";
  std::ofstream tor_500sample(filename);

  std::vector<Point_3> samplePoints = n_torus_sample_points(500, 1, 3, 3109874);

  for (Point_3 p: samplePoints) {
    tor_500sample << p.x() << ", " << p.y() << ", " << p.z() << std::endl;
  }

  
  return 0; 
}
