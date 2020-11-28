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
    Vecteur normalize();

    Vecteur operator -() {
        return Vecteur(-x, -y, -z);
    }

    Vecteur operator *(float fact) {
        return Vecteur(x*fact, y*fact, z*fact);
    }
};

#endif // VECTEUR_H
