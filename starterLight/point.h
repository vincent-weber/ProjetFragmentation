#ifndef POINT_H
#define POINT_H


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

#endif // POINT_H
