#ifndef TETRAEDRE_H
#define TETRAEDRE_H

#include "point.h"
#include <array>

class Tetraedre
{
public:
    std::array<Point*, 4> points;
    std::array<Tetraedre*, 4> tetra_adj;

    //Pour respecter critere boule vide
    Point circumsphere_center;

    //Pour l'algorithme walk il le faut apparemment
    Point insphere_center;

    Tetraedre(Point* p1, Point* p2, Point* p3, Point* p4);
    Tetraedre();
};

#endif // TETRAEDRE_H
