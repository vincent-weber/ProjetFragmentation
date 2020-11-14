#ifndef TETRAEDRE_H
#define TETRAEDRE_H

#include "point.h"
#include "vecteur.h"
#include <array>

class Tetraedre
{
public:
    std::array<Point*, 4> points;
    std::array<Tetraedre*, 4> tetra_adj;
    //normales[0] = normale a la face formée par les points points[0], [1] et [2]
    //normales[1] = normale a la face formée par les points points[1], [2] et [3]
    //normales[2] = normale a la face formée par les points points[2], [3] et [0]
    //normales[3] = normale a la face formée par les points points[3], [0] et [1]
    std::array<Vecteur, 4> normales = {};

    //Pour respecter critere boule vide
    Point circumsphere_center;

    //Pour l'algorithme walk il le faut apparemment
    Point insphere_center;

    Tetraedre(Point* p1, Point* p2, Point* p3, Point* p4);
    Tetraedre();
};

#endif // TETRAEDRE_H
