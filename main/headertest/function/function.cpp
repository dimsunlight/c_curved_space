#include <iostream> 
#include <math.h>

float add(float f1, float f2) {
	return f1 + f2;
}

float subtract(float f1, float f2) {
	return f1 - f2;
}

float mumboJumbo (float f1, float f2) {
	float weirdStuff = add(f1,f2)*subtract(f1,f2) - add(subtract(f1,f2),subtract(f1,f2));
	return weirdStuff;
}
