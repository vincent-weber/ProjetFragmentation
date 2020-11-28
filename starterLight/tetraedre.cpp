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
            orient_faces[i] = false;
        } else {
            normales[i] = normale;
            orient_faces[i] = true;
        }
    }

    float norm1 = normales[0].norm();
    float norm2 = normales[1].norm();
    float norm3 = normales[2].norm();
    float norm4 = normales[3].norm();

    insphere_center = (*points[3] * norm1 + *points[0] * norm2 + *points[1] * norm3 + *points[2] * norm4) / (norm1 + norm2 + norm3 + norm4);
}

Tetraedre::Tetraedre() {}

bool Tetraedre::contientPoint(Point* p) {
    return (p == points[0] || p == points[1] || p == points[2] || p == points[3]);
}

bool Tetraedre::isPointInSphere(Point* p){
    std::vector<float> a;
    a.push_back(this->points[0]->x);
    a.push_back(this->points[0]->y);
    a.push_back(this->points[0]->z);
    a.push_back(1);

    std::vector<float> b;
    b.push_back(this->points[1]->x);
    b.push_back(this->points[1]->y);
    b.push_back(this->points[1]->z);
    b.push_back(1);


    std::vector<float> c;
    c.push_back(this->points[2]->x);
    c.push_back(this->points[2]->y);
    c.push_back(this->points[2]->z);
    c.push_back(1);

    std::vector<float> d;
    d.push_back(this->points[3]->x);
    d.push_back(this->points[3]->y);
    d.push_back(this->points[3]->z);
    d.push_back(1);

    std::vector<std::vector<float>> orientMatrix;
    orientMatrix.push_back(a);
    orientMatrix.push_back(b);
    orientMatrix.push_back(c);
    orientMatrix.push_back(d);
    float detOrient = calcDeterminant(orientMatrix);

    a.at(3) = pow(this->points[0]->x, 2) + pow(this->points[0]->y, 2) + pow(this->points[0]->z, 2);
    a.push_back(1);

    b.at(3) = pow(this->points[1]->x, 2) + pow(this->points[1]->y, 2) + pow(this->points[1]->z, 2);
    b.push_back(1);

    c.at(3) = pow(this->points[2]->x, 2) + pow(this->points[2]->y, 2) + pow(this->points[2]->z, 2);
    c.push_back(1);

    d.push_back(pow(this->points[3]->x, 2) + pow(this->points[3]->y, 2) + pow(this->points[3]->z, 2));
    d.push_back(1);

    std::vector<float> mp;
    mp.push_back(p->x);
    mp.push_back(p->y);
    mp.push_back(p->z);
    mp.push_back(pow(p->x, 2) + pow(p->y, 2) + pow(p->z, 2));
    mp.push_back(1);

    std::vector<std::vector<float>> inSphereMatrix;
    if(detOrient > 0){
        inSphereMatrix.push_back(a);
        inSphereMatrix.push_back(b);
    }
    else{
        inSphereMatrix.push_back(b);
        inSphereMatrix.push_back(a);
    }
    inSphereMatrix.push_back(c);
    inSphereMatrix.push_back(d);
    inSphereMatrix.push_back(mp);

    float detInSphere = calcDeterminant(inSphereMatrix);
    qDebug() << "DETERMINANT IN SPHERE : " <<  detInSphere;
    if(detInSphere > 0)
        return true;
    else
        return false;
}


QDebug operator <<(QDebug os, const Tetraedre& t) {
    for (int i = 0 ; i < 4 ; ++i) {
        os << *t.points[i] << '\n';
    }
    return os;
}
