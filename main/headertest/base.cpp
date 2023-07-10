#include <iostream>
#include "function.h"

int main() {
  float f1 = 23.9;
  float f2 = 5.7;  
  std::cout << "The floats are " << f1 << " and " << f2 << std::endl;
  float completeMess = mumboJumbo(f1,f2);
  std::cout << "the post-function quantity is " << completeMess << std::endl;

  return 0;
}
