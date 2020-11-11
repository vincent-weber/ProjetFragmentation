#ifndef VECTEUR_H
#define VECTEUR_H

#include "point.h"
#include <cmath>

class Vecteur
{
public:
    float x;
    float y;
    float z;

    Vecteur(Point p1, Point p2);
    Vecteur(float x, float y, float z);
    Vecteur();
    float norm();
    float dot_product(Vecteur v);
    Vecteur cross_product(Vecteur v);

    Vecteur operator -() {
        return Vecteur(-x, -y, -z);
    }
};

#endif // VECTEUR_H
