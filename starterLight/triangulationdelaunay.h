#ifndef TRIANGULATIONDELAUNAY_H
#define TRIANGULATIONDELAUNAY_H

#include "float.h"
#include <stack>
#include <vector>
#include <forward_list>
#include "point.h"
#include "tetraedre.h"

class TriangulationDelaunay
{
public:

    std::stack<Point*> remaining_points;
    std::forward_list<Tetraedre*> tetraedres;
    //scale_factor au moins 50 je dirais
    float big_tetra_scale_factor = 20.0f;

    //est-ce que l'algo return une liste de points (suite de 4 points = 1 tetraedre) ou de tetraedres ?
    TriangulationDelaunay(std::vector<Point>* tetra_ret, std::vector<Point> points);

    void triangulation();

    Tetraedre* tetra_containing_point_walk(Point* point);

    //Flips standards
    std::array<Tetraedre*, 2> flip32(Tetraedre* t1, Tetraedre* t2, Tetraedre* t3);
    std::array<Tetraedre*, 3> flip23(Tetraedre* t1, Tetraedre* t2);
    std::array<Tetraedre*, 4> flip14(Tetraedre* t1);

    //Flips pour cas degeneres
    std::array<Tetraedre*, 4> flip44(Tetraedre* t1, Tetraedre* t2, Tetraedre* t3, Tetraedre* t4);
};

#endif // TRIANGULATIONDELAUNAY_H
