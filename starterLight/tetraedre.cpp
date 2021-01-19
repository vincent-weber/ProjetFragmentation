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
    computeCircumCenter();
}

Tetraedre::Tetraedre() {}

void Tetraedre::computeCircumCenter() {
    std::vector<float> aO;
    aO.push_back(this->points[0]->x);
    aO.push_back(this->points[0]->y);
    aO.push_back(this->points[0]->z);
    aO.push_back(1);

    std::vector<float> bO;
    bO.push_back(this->points[1]->x);
    bO.push_back(this->points[1]->y);
    bO.push_back(this->points[1]->z);
    bO.push_back(1);

    std::vector<float> cO;
    cO.push_back(this->points[2]->x);
    cO.push_back(this->points[2]->y);
    cO.push_back(this->points[2]->z);
    cO.push_back(1);

    std::vector<float> dO;
    dO.push_back(this->points[3]->x);
    dO.push_back(this->points[3]->y);
    dO.push_back(this->points[3]->z);
    dO.push_back(1);

    std::vector<std::vector<float>> alphaMatrix;
    alphaMatrix.push_back(aO);
    alphaMatrix.push_back(bO);
    alphaMatrix.push_back(cO);
    alphaMatrix.push_back(dO);
    float alpha = calcDeterminant(alphaMatrix);

    std::vector<std::vector<float>> DX;
    std::vector<std::vector<float>> DY;
    std::vector<std::vector<float>> DZ;

    float a_squared = points[0]->x * points[0]->x + points[0]->y * points[0]->y + points[0]->z * points[0]->z;
    float b_squared = points[1]->x * points[1]->x + points[1]->y * points[1]->y + points[1]->z * points[1]->z;
    float c_squared = points[2]->x * points[2]->x + points[2]->y * points[2]->y + points[2]->z * points[2]->z;
    float d_squared = points[3]->x * points[3]->x + points[3]->y * points[3]->y + points[3]->z * points[3]->z;

    std::vector<float> DX1;
    DX1.push_back(a_squared); DX1.push_back(points[0]->y); DX1.push_back(points[0]->z); DX1.push_back(1);
    std::vector<float> DX2;
    DX2.push_back(b_squared); DX2.push_back(points[1]->y); DX2.push_back(points[1]->z); DX2.push_back(1);
    std::vector<float> DX3;
    DX3.push_back(c_squared); DX3.push_back(points[2]->y); DX3.push_back(points[2]->z); DX3.push_back(1);
    std::vector<float> DX4;
    DX4.push_back(d_squared); DX4.push_back(points[3]->y); DX4.push_back(points[3]->z); DX4.push_back(1);
    DX.push_back(DX1); DX.push_back(DX2); DX.push_back(DX3); DX.push_back(DX4);

    std::vector<float> DY1;
    DY1.push_back(a_squared); DY1.push_back(points[0]->x); DY1.push_back(points[0]->z); DY1.push_back(1);
    std::vector<float> DY2;
    DY2.push_back(b_squared); DY2.push_back(points[1]->x); DY2.push_back(points[1]->z); DY2.push_back(1);
    std::vector<float> DY3;
    DY3.push_back(c_squared); DY3.push_back(points[2]->x); DY3.push_back(points[2]->z); DY3.push_back(1);
    std::vector<float> DY4;
    DY4.push_back(d_squared); DY4.push_back(points[3]->x); DY4.push_back(points[3]->z); DY4.push_back(1);
    DY.push_back(DY1); DY.push_back(DY2); DY.push_back(DY3); DY.push_back(DY4);

    std::vector<float> DZ1;
    DZ1.push_back(a_squared); DZ1.push_back(points[0]->x); DZ1.push_back(points[0]->y); DZ1.push_back(1);
    std::vector<float> DZ2;
    DZ2.push_back(b_squared); DZ2.push_back(points[1]->x); DZ2.push_back(points[1]->y); DZ2.push_back(1);
    std::vector<float> DZ3;
    DZ3.push_back(c_squared); DZ3.push_back(points[2]->x); DZ3.push_back(points[2]->y); DZ3.push_back(1);
    std::vector<float> DZ4;
    DZ4.push_back(d_squared); DZ4.push_back(points[3]->x); DZ4.push_back(points[3]->y); DZ4.push_back(1);
    DZ.push_back(DZ1); DZ.push_back(DZ2); DZ.push_back(DZ3); DZ.push_back(DZ4);

    float DetDX = calcDeterminant(DX);
    float DetDY = -calcDeterminant(DY);
    float DetDZ = calcDeterminant(DZ);
    circumsphere_center = Point(DetDX, DetDY, DetDZ) / (2*alpha);
    circumsphere_radius = sqrt(pow(circumsphere_center.x - points[0]->x, 2) + pow(circumsphere_center.y - points[0]->y, 2) + pow(circumsphere_center.z - points[0]->z, 2));

}

bool Tetraedre::contientPoint(Point* p) {
    return (p == points[0] || p == points[1] || p == points[2] || p == points[3]);
}

//DORENAVANT INUTILISE AU PROFIT DU CALCUL DU CENTRE DU CERCLE CIRONSCRIT DIRECT
bool Tetraedre::isPointInSphere(Point* p){
    std::vector<float> aO;
    aO.push_back(this->points[0]->x);
    aO.push_back(this->points[0]->y);
    aO.push_back(this->points[0]->z);
    aO.push_back(1);

    std::vector<float> bO;
    bO.push_back(this->points[1]->x);
    bO.push_back(this->points[1]->y);
    bO.push_back(this->points[1]->z);
    bO.push_back(1);

    std::vector<float> cO;
    cO.push_back(this->points[2]->x);
    cO.push_back(this->points[2]->y);
    cO.push_back(this->points[2]->z);
    cO.push_back(1);

    std::vector<float> dO;
    dO.push_back(this->points[3]->x);
    dO.push_back(this->points[3]->y);
    dO.push_back(this->points[3]->z);
    dO.push_back(1);

    std::vector<float> mpO;
    mpO.push_back(p->x);
    mpO.push_back(p->y);
    mpO.push_back(p->z);
    mpO.push_back(1);

    std::vector<std::vector<float>> orientMatrix;
    /*Vecteur v1(*this->points[0],*this->points[1]);
    Vecteur v2(*this->points[0],*this->points[2]);
    Vecteur vn = v1.cross_product((v2));*/
    if (orient_faces[0]) {
        orientMatrix.push_back(aO);
        orientMatrix.push_back(bO);
    }
    else {
        orientMatrix.push_back(bO);
        orientMatrix.push_back(aO);
    }
    orientMatrix.push_back(cO);
    orientMatrix.push_back(dO);
    float detOrient = calcDeterminant(orientMatrix);
    qDebug() << "DET ORIENT : " << detOrient;

    std::vector<float> a = aO;
    a.at(3) = points[0]->x * points[0]->x + points[0]->y * points[0]->y + points[0]->z * points[0]->z;
    a.push_back(1);

    std::vector<float> b = bO;
    b.at(3) = points[1]->x * points[1]->x + points[1]->y * points[1]->y + points[1]->z * points[1]->z;
    b.push_back(1);

    std::vector<float> c = cO;
    c.at(3) = points[2]->x * points[2]->x + points[2]->y * points[2]->y + points[2]->z * points[2]->z;
    c.push_back(1);

    std::vector<float> d = dO;
    d.at(3) = points[3]->x * points[3]->x + points[3]->y * points[3]->y + points[3]->z * points[3]->z;
    d.push_back(1);

    std::vector<float> mp = mpO;
    mp.at(3) = p->x * p->x + p->y * p->y + p->z * p->z;
    mp.push_back(1);

    std::vector<std::vector<float>> inSphereMatrix;
    std::vector<std::vector<float>> orientMatrixDebug;
    if(detOrient > 0){
        inSphereMatrix.push_back(a);
        inSphereMatrix.push_back(b);
        orientMatrixDebug.push_back(aO);
        orientMatrixDebug.push_back(bO);
    }
    else{
        inSphereMatrix.push_back(b);
        inSphereMatrix.push_back(a);
        orientMatrixDebug.push_back(bO);
        orientMatrixDebug.push_back(aO);
    }
    inSphereMatrix.push_back(c);
    inSphereMatrix.push_back(d);
    inSphereMatrix.push_back(mp);
    orientMatrixDebug.push_back(cO);
    orientMatrixDebug.push_back(dO);

    float detInSphere = calcDeterminant(inSphereMatrix);
    if(detInSphere > 0)
        return true;
    else
        return false;
}

bool Tetraedre::write_tetra_to_file(const std::string& filename, Tetraedre* tetra) {

}



QDebug operator <<(QDebug os, const Tetraedre& t) {
    os << "Tetra courant : " << &t << '\n';
    for (int i = 0 ; i < 4 ; ++i) {
        os << *t.points[i] << '\n';
    }
    os << "Tetra adjacents : " << '\n';
    for (int i = 0 ; i < 4 ; ++i) {
        os << t.tetra_adj[i] << '\n';
    }
    return os;
}
