#include "vecteur.h"

Vecteur::Vecteur(Point p1, Point p2)
{
    this->x = p2.x - p1.x;
    this->y = p2.y - p1.y;
    this->z = p2.z - p1.z;
}

Vecteur::Vecteur(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

Vecteur::Vecteur() {}

float Vecteur::norm() {
    return sqrt(x*x + y*y + z*z);
}

float Vecteur::dot_product(Vecteur v) {
    return this->x * v.x + this->y * v.y + this->z * v.z;
}

Vecteur Vecteur::cross_product(Vecteur v) {
    float x = this->y * v.z - this->z * v.y;
    float y = this->z * v.x - this->x * v.z;
    float z = this->x * v.y - this->y * v.x;
    return Vecteur(x,y,z);
}
