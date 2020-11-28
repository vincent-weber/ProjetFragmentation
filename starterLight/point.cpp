#include "point.h"

Point::Point(float x, float y, float z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

Point::Point() {}

QDebug operator <<(QDebug os, const Point& p) {
    os.nospace() << p.x << " - " << p.y << " - " << p.z;
    return os;
}
