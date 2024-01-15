/*--------------------------------------------------
  Gprof Example 2
  --------------------------------------------------
  main
    |-- func1 -- func3
    |-- func2 
  --------------------------------------------------
*/
#include <cstdlib>
#include <stdio.h>
using namespace std;
//declaration of functions
int func1();
int func2();
int func3();

int func1(void) {
    int i=0,g=0;
    while(i++<10000) {
        g+=i;
        func3();
    }
    return g;
}

int func2(void) {
    int i=0,g=0;
    while(i++<40000) {
        g+=i;
    }
    return g;
}


int func3() {
    int i=0,g=0;
    while(i++ < 1000) {
        g++;
    }
    return g;
}


int main(int argc, char** argv) {
    int iterations = 1000;
    printf("Number of iterations = %d\n", iterations);
    while(iterations--) {
        printf("the iteration No. %d\n", 1000-iterations);
        func1();
        func2();
    }
}  
