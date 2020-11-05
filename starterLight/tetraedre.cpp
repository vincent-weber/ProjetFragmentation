#include "tetraedre.h"

Tetraedre::Tetraedre(Point* p1, Point* p2, Point* p3, Point* p4)
{
    this->points[0] = p1;
    this->points[1] = p2;
    this->points[2] = p3;
    this->points[3] = p4;

    tetra_adj[0] = nullptr;
    tetra_adj[1] = nullptr;
    tetra_adj[2] = nullptr;
    tetra_adj[3] = nullptr;
}

Tetraedre::Tetraedre() {}
