#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include "math.h"
#include "vecteur.h"
#include "point.h"

int calcDeterminant(std::vector<std::vector<float>> Matrix);
bool intersect_droite_triangle(Vecteur& vec_direct_droite, Point& p_droite, Point& p1, Point& p2, Point& p3);


#endif // UTIL_H
