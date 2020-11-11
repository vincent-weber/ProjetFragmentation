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

    for (int i = 0 ; i < 4 ; ++i) {
        Vecteur vf1(*points[i], *points[(i+1)%4]);
        Vecteur vf2(*points[i], *points[(i+2)%4]);
        Vecteur v3(*points[i], *points[(i+3)%4]);

        Vecteur normale = vf1.cross_product(vf2);
        float sign = normale.dot_product(v3);

        if (sign > 0) {
            normales[i] = -normale;
        } else {
            normales[i] = normale;
        }
    }

    float norm1 = normales[0].norm();
    float norm2 = normales[1].norm();
    float norm3 = normales[2].norm();
    float norm4 = normales[3].norm();

    insphere_center = (*points[3] * norm1 + *points[0] * norm2 + *points[1] * norm3 + *points[2] * norm4) / (norm1 + norm2 + norm3 + norm4);
}

Tetraedre::Tetraedre() {}
