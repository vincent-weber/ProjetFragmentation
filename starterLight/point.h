#ifndef POINT_H
#define POINT_H

#include <QDebug>

class Point
{
public:
    float x;
    float y;
    float z;

    Point(float x, float y, float z);
    Point();

    Point operator *(float fact) {
        return Point(x*fact, y*fact, z*fact);
    }

    Point operator +(Point p) {
        return Point(x + p.x, y + p.y, z + p.z);
    }

    Point operator /(float div) {
        return Point(x / div, y / div, z / div);
    }
};

QDebug operator <<(QDebug os, const Point& p);

#endif // POINT_H
