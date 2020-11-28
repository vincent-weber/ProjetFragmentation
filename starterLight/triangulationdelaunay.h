#ifndef TRIANGULATIONDELAUNAY_H
#define TRIANGULATIONDELAUNAY_H

#include "float.h"
#include <stack>
#include <vector>
#include <forward_list>
#include <unordered_set>
#include "point.h"
#include "tetraedre.h"
#include "util.h"
#include <QDebug>

class TriangulationDelaunay
{
public:

    unsigned flip14count = 0;
    unsigned flip23count = 0;
    unsigned flip32count = 0;

    std::stack<Point*> remaining_points;
    std::forward_list<Tetraedre*> tetraedres;
    //scale_factor au moins 50 je dirais
    float big_tetra_scale_factor = 28.0f;

    std::array<Point*,4> p_tetra_englob;
    unsigned tetra_count = 0;

    //est-ce que l'algo return une liste de points (suite de 4 points = 1 tetraedre) ou de tetraedres ?
    TriangulationDelaunay(std::vector<Point*>& points);

    void triangulation();

    void remplacer_tetraedre_adjacent(Tetraedre* tetra_a_modif_adj, Tetraedre* tetra_adj_a_remplacer, Tetraedre* tetra_adj_remplacant);
    Tetraedre* tetraedre_adjacent_face(Tetraedre* t, Point* p1, Point* p2, Point* p3);
    Tetraedre* tetra_containing_point_walk(Point* point);

    //Flips standards
    std::array<Tetraedre*, 2> flip32(std::array<Tetraedre*, 3>);
    std::array<Tetraedre*, 3> flip23(std::array<Tetraedre*, 2>);
    std::array<Tetraedre*, 4> flip14(Tetraedre* t1, Point* p);

    //Flips pour cas degeneres
    std::array<Tetraedre*, 4> flip44(std::array<Tetraedre*, 4>);
};

#endif // TRIANGULATIONDELAUNAY_H
