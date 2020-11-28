#include "triangulationdelaunay.h"

TriangulationDelaunay::TriangulationDelaunay(std::vector<Point*>& points)
{

    float min_x = FLT_MAX, min_y = FLT_MAX, min_z = FLT_MAX;
    float max_x = FLT_MIN, max_y = FLT_MIN, max_z = FLT_MIN;
    for (Point* point : points) {
        remaining_points.push(point);
        if (point->x < min_x) min_x = point->x;
        else if (point->x > max_x) max_x = point->x;
        if (point->y < min_y) min_y = point->y;
        else if (point->y > max_y) max_y = point->y;
        if (point->z < min_z) min_z = point->z;
        else if (point->z > max_z) max_z = point->z;
    }

    /*for (int i = 0 ; i < points.size() ; ++i) {
        qDebug() << points[i];
    }*/

    float center_x = (max_x + min_x) / 2;
    float dist_center_x = (max_x - min_x) / 2;
    float x_box_min = center_x - (big_tetra_scale_factor * dist_center_x);
    float x_box_max = -x_box_min;

    float center_y = (max_y + min_y) / 2;
    float dist_center_y = (max_y - min_y) / 2;
    float y_box_min = center_y - (big_tetra_scale_factor * dist_center_y);
    float y_box_max = -y_box_min;

    float center_z = (max_z + min_z) / 2;
    float dist_center_z = (max_z - min_z) / 2;
    float z_box_min = center_z - (big_tetra_scale_factor * dist_center_z);
    float z_box_max = -z_box_min;


    Point* p1 = new Point(center_x, y_box_max, center_z);
    Point* p2 = new Point(x_box_min, y_box_min, z_box_min);
    Point* p3 = new Point(x_box_max, y_box_min, z_box_min);
    Point* p4 = new Point(center_x, y_box_min, z_box_max);

    //qDebug() << p1 << p2 << p3 << p4;

    Tetraedre* bigT = new Tetraedre(p1, p2, p3, p4);
    p_tetra_englob[0] = p1;
    p_tetra_englob[1] = p2;
    p_tetra_englob[2] = p3;
    p_tetra_englob[3] = p4;

    std::vector<Point> points_ret;

    tetraedres.push_front(bigT);
    tetra_count = 1;

    triangulation();

}

void TriangulationDelaunay::remplacer_tetraedre_adjacent(Tetraedre* tetra_a_modif_adj, Tetraedre* tetra_adj_a_remplacer, Tetraedre* tetra_adj_remplacant) {
    if (tetra_a_modif_adj ==  nullptr) return;
    for (int i = 0 ; i < 4 ; ++i) {
        if (tetra_a_modif_adj->tetra_adj[i] == tetra_adj_a_remplacer) {
            tetra_a_modif_adj->tetra_adj[i] = tetra_adj_remplacant;
            return;
        }
    }
}

Tetraedre* TriangulationDelaunay::tetraedre_adjacent_face(Tetraedre* t, Point* p1, Point* p2, Point* p3) {
    if (t == nullptr || p1 == nullptr || p2 == nullptr || p3 == nullptr) return nullptr;
    for (Tetraedre* t_adj : t->tetra_adj) {
        if (t_adj != nullptr) {
            bool ok1 = false, ok2 = false, ok3 = false;
            for (Point* p_t_adj : t_adj->points) {
                if (p_t_adj == p1) {
                    ok1 = true;
                }
                if (p_t_adj == p2) {
                    ok2 = true;
                }
                if (p_t_adj == p3) {
                    ok3 = true;
                }
            }
            if (ok1 && ok2 && ok3) {
                return t_adj;
            }
            ok1 = false;
            ok2 = false;
            ok3 = false;
        }
    }
    return nullptr;
}

std::array<Tetraedre*, 4> TriangulationDelaunay::flip14(Tetraedre *old_tetra, Point* p){
    std::array<Tetraedre*, 4> newTetras;
    std::array<Tetraedre*, 4> adj;
    for(int i = 0; i < 4; i++) {
        Tetraedre* nouvTetra = new Tetraedre(p, old_tetra->points[i%4], old_tetra->points[(i+1)%4], old_tetra->points[(i+2)%4]);
        newTetras[i] = nouvTetra;
        adj[i] = tetraedre_adjacent_face(old_tetra, old_tetra->points[i%4], old_tetra->points[(i+1)%4], old_tetra->points[(i+2)%4]);
    }
    for(int i = 0; i < 4; i++){
        newTetras[i]->tetra_adj = {newTetras[(i+1)%4], newTetras[(i+2)%4], newTetras[(i+3)%4], adj[i]};
    }
    //TODO : A priori, on peut fusionner cette boucle for avec la précédente, a tester
    for(int i = 0; i < 4; i++){
        remplacer_tetraedre_adjacent(adj[i], old_tetra, newTetras[i]);
        tetraedres.push_front(newTetras[i]);
    }
    tetraedres.remove(old_tetra);
    delete old_tetra;

    tetra_count += 3;
    return newTetras;
}

std::array<Tetraedre*, 3> TriangulationDelaunay::flip23(std::array<Tetraedre*, 2> old_tetras) {
    //on sépare les 5 sommets différents en 2 catégories.
    Point* diff_point1 = nullptr;
    Point* diff_point2 = nullptr;
    std::array<Point*, 3> common_points;
    std::map<Point*, int> nb_occurrences_sommets;
    for (Tetraedre* tetra : old_tetras) {
        for (int i = 0 ; i < 4 ; ++i) {
            nb_occurrences_sommets.insert(std::pair<Point*, int>(tetra->points[i], 0));
            ++nb_occurrences_sommets.at(tetra->points[i]);
        }
    }

    int i = 0;
    for (auto it = nb_occurrences_sommets.begin() ; it != nb_occurrences_sommets.end() ; ++it) {
        if (it->second == 1 && diff_point1 == nullptr) {
            diff_point1 = it->first;
        }
        else if (it->second == 1 && diff_point2 == nullptr) {
            diff_point2 = it->first;
        }
        else {
            common_points[i++] = it->first;
        }
    }

    Tetraedre* tetra1 = new Tetraedre(diff_point1, diff_point2, common_points[0], common_points[1]);
    Tetraedre* tetra2 = new Tetraedre(diff_point1, diff_point2, common_points[1], common_points[2]);
    Tetraedre* tetra3 = new Tetraedre(diff_point1, diff_point2, common_points[2], common_points[0]);
    std::array<Tetraedre*, 3> tetras_ret = {tetra1, tetra2, tetra3};

    //parcours des 3 nouveaux tetraedres : i = indice nouveau tetraedre
    for (int i = 0 ; i < 3 ; ++i) {
        Tetraedre* tetra_courant = tetras_ret[i];
        std::vector<Tetraedre*> nouv_tetras_adj;
        //parcours des 2 faces du tetra courant qui nous intéressent (faces extérieures) : j = indice face
        for (int j = 0 ; j < 2 ; ++j) {
            Point* p1_face = tetra_courant->points[j];
            Point* p2_face = tetra_courant->points[2];
            Point* p3_face = tetra_courant->points[3];
            //parcours des 2 tetraedres anciens : k = indice ancien tetraedre
            //On cherche l'ancien tetraedre qui contient la face en cours
            for (int k = 0 ; k < 2 ; ++k) {
                if (old_tetras[k]->contientPoint(p1_face) && old_tetras[k]->contientPoint(p2_face) && old_tetras[k]->contientPoint(p3_face)) {
                    //update des adjacences
                    Tetraedre* tetra_adj_a_modif = tetraedre_adjacent_face(old_tetras[k], p1_face, p2_face, p3_face);
                    nouv_tetras_adj.push_back(tetra_adj_a_modif);
                    remplacer_tetraedre_adjacent(tetra_adj_a_modif/*tetra sur lequel remplacer*/, old_tetras[k]/*remplacé*/, tetra_courant/*remplacant*/);
                    break;
                }
            }
        }
        tetras_ret[i]->tetra_adj = {tetras_ret[(i+1)%3], tetras_ret[(i+2)%3], nouv_tetras_adj[0], nouv_tetras_adj[1]};
    }

    tetraedres.remove(old_tetras[0]); tetraedres.remove(old_tetras[1]);
    delete old_tetras[0]; delete old_tetras[1];
    tetraedres.push_front(tetras_ret[0]); tetraedres.push_front(tetras_ret[1]); tetraedres.push_front(tetras_ret[2]);

    tetra_count += 1;
    return tetras_ret;
}

std::array<Tetraedre*, 2> TriangulationDelaunay::flip32(std::array<Tetraedre*, 3> old_tetras) {
    //on sépare les 5 sommets différents en 2 catégories.
    Point* diff_point1 = nullptr;
    Point* diff_point2 = nullptr;
    std::array<Point*, 3> common_points;
    std::map<Point*, int> nb_occurrences_sommets;
    for (Tetraedre* tetra : old_tetras) {
        for (int i = 0 ; i < 4 ; ++i) {
            nb_occurrences_sommets.insert(std::pair<Point*, int>(tetra->points[i], 0));
            ++nb_occurrences_sommets.at(tetra->points[i]);
        }
    }

    int i = 0;
    for (auto it = nb_occurrences_sommets.begin() ; it != nb_occurrences_sommets.end() ; ++it) {
        if (it->second == 3 && diff_point1 == nullptr) {
            diff_point1 = it->first;
        }
        else if (it->second == 3 && diff_point2 == nullptr) {
            diff_point2 = it->first;
        }
        else {
            common_points[i++] = it->first;
        }
    }

    Tetraedre* tetra1 = new Tetraedre(common_points[0], common_points[1], common_points[2], diff_point1);
    Tetraedre* tetra2 = new Tetraedre(common_points[0], common_points[1], common_points[2], diff_point2);
    std::array<Tetraedre*, 2> tetras_ret = {tetra1, tetra2};

    //parcours des 2 nouveaux tetraedres : i = indice nouveau tetraedre
    for (int i = 0 ; i < 2 ; ++i) {
        Tetraedre* tetra_courant = tetras_ret[i];
        std::vector<Tetraedre*> nouv_tetras_adj;
        //parcours des 3 faces du tetra courant qui nous intéressent (faces extérieures, pas la face commune) : j = indice face
        for (int j = 0 ; j < 3 ; ++j) {
            Point* p1_face = tetra_courant->points[3]; //point_diff
            Point* p2_face = tetra_courant->points[j];
            Point* p3_face = tetra_courant->points[(j+1)%3];
            //parcours des 3 tetraedres anciens : k = indince ancien tetraedre
            //On cherche l'ancien tetraedre qui contient la face en cours
            for (int k = 0 ; k < 3 ; ++k) {
                if (old_tetras[k]->contientPoint(p1_face) && old_tetras[k]->contientPoint(p2_face) && old_tetras[k]->contientPoint(p3_face)) {
                    //update des adjacences
                    Tetraedre* tetra_adj_a_modif = tetraedre_adjacent_face(old_tetras[k], p1_face, p2_face, p3_face);
                    nouv_tetras_adj.push_back(tetra_adj_a_modif);
                    remplacer_tetraedre_adjacent(tetra_adj_a_modif/*tetra sur lequel remplacer*/, old_tetras[k]/*remplacé*/, tetra_courant/*remplacant*/);
                    break;
                }
            }
        }
        tetras_ret[i]->tetra_adj = {tetras_ret[(i+1)%2], nouv_tetras_adj[0], nouv_tetras_adj[1], nouv_tetras_adj[2]};
    }

    tetraedres.remove(old_tetras[0]); tetraedres.remove(old_tetras[1]); tetraedres.remove(old_tetras[2]);
    delete old_tetras[0]; delete old_tetras[1]; delete old_tetras[2];
    tetraedres.push_front(tetras_ret[0]); tetraedres.push_front(tetras_ret[1]);

    tetra_count -= 1;
    return tetras_ret;
}

Tetraedre* TriangulationDelaunay::tetra_containing_point_walk(Point* point) {
    //algorithme walk
    //TODO : boucles infinies = moche, voir ce qu'on fait
    //TODO : voir ce qu'on fait pour les tetras deja visites
    std::vector<Tetraedre*> visited_tetras;
    Tetraedre* tetra = tetraedres.front();
    while(true) {
        visited_tetras.push_back(tetra);
        //check si le point est dans le tetra courant
        //etant donné que tetra->normales[1] est TOUJOURS la normale à la face qui ne contient pas le point tetra->points[0],
        //le dot product entre v_test et normales[1] sera toujours positif. Alors on créé un v_test2 issu d'un des 3 autres sommets
        //autre idée si ca marche pas : on fait le dot product entre la normale des faces et le vecteur issu du sommet opposé
        //aux faces qui va vers le insphere center, et si les 4 dot products sont positifs, on return le tetra.
        Vecteur v_test(*tetra->points[0], tetra->insphere_center);
        Vecteur v_test2(*tetra->points[1], tetra->insphere_center);
        float dot1 = tetra->normales[0].dot_product(v_test);
        float dot2 = tetra->normales[1].dot_product(v_test2);
        float dot3 = tetra->normales[2].dot_product(v_test);
        float dot4 = tetra->normales[3].dot_product(v_test);
        if (dot1 < 0 && dot2 < 0 && dot3 < 0 && dot4 < 0) {
            return tetra;
        }

        //le tetraedre courant n'est pas le bon, on cherche la face vers laquelle se diriger ensuite
        Vecteur v(tetra->insphere_center, *point);
        float max_dot = FLT_MIN;
        int ind_face = -1;
        for (int i = 0 ; i < 4 ; ++i) {
            //je suis quasi sur qu'il faut normaliser la normale
            float dot = tetra->normales[i].normalize().dot_product(v);
            if (dot > max_dot) {
                max_dot = dot;
                ind_face = i;
            }
        }

        //on memorise les points de la face pour laquelle il faut trouver le tetra adjacent
        Point* p_face_1 = tetra->points[ind_face];
        Point* p_face_2 = tetra->points[(ind_face+1)%4];
        Point* p_face_3 = tetra->points[(ind_face+2)%4];

        //on cherche le tetra adjacent par la face trouvee precedemment
        for (int ind_tetra_adj = 0 ; ind_tetra_adj < 4 ; ++ind_tetra_adj) {
            if (tetra->tetra_adj[ind_tetra_adj] == nullptr) {
                continue;
            }
            Tetraedre* tetra_adj = tetra->tetra_adj[ind_tetra_adj];
            int common_points = 0;
            if (tetra_adj->contientPoint(p_face_1)) ++common_points;
            if (tetra_adj->contientPoint(p_face_2)) ++common_points;
            if (tetra_adj->contientPoint(p_face_3)) ++common_points;
            //si les 3 points de la face appartiennent au triangle adjacent courant, c'est le bon
            if (common_points == 3) {
                tetra = tetra_adj;
            }
        }
    }


    return nullptr;
}

Tetraedre* tetra_adj_2_points(Tetraedre* t1, Point* p1, Tetraedre* t2, Point* p2) {
    std::vector<Tetraedre*> adj_t1;
    std::vector<Tetraedre*> adj_t2;
    for (int i = 0 ; i < 4 ; ++i) {
        if (t1->tetra_adj[i] != nullptr && t1->tetra_adj[i]->contientPoint(p1)) {
            adj_t1.push_back(t1->tetra_adj[i]);
        }
        if (t2->tetra_adj[i] != nullptr && t2->tetra_adj[i]->contientPoint(p2)) {
            adj_t2.push_back(t2->tetra_adj[i]);
        }
    }

    for (int i = 0 ; i < 3 ; ++i) {
        for (int j = 0 ; j < 3 ; ++j) {
            if (adj_t1[i] == adj_t2[j]) {
                return adj_t1[i];
            }
        }
    }
    return nullptr;
}

void TriangulationDelaunay::triangulation() {
    while (!remaining_points.empty()) {
        Point* p = remaining_points.top();
        remaining_points.pop();
        Tetraedre* current_tetra = tetra_containing_point_walk(p);
        std::array<Tetraedre*, 4> tetras_new = flip14(current_tetra, p);
        ++flip14count;
        std::stack<Tetraedre*> tetras_to_check;
        for (int i = 0 ; i < 4 ; ++i) {
            tetras_to_check.push(tetras_new[i]);
        }

        while (!tetras_to_check.empty()) {
            Tetraedre* t = tetras_to_check.top();
            tetras_to_check.pop();
            std::array<Point*, 3> abc = {nullptr, nullptr, nullptr};
            for (int i = 0 ; i < 4 ; ++i) {
                if (t->points[i] == p) {
                    abc[0] = t->points[(i+1)%4];
                    abc[1] = t->points[(i+2)%4];
                    abc[2] = t->points[(i+3)%4];
                    break;
                }
            }
            if (abc[0] == nullptr) {
                qDebug() << "Erreur : le tetraedre courant ne contient pas le point p.";
                //exit(-1);
            }

            Tetraedre* ta = tetraedre_adjacent_face(t, abc[0], abc[1], abc[2]);
            if (ta == nullptr) continue;
            Point* d = nullptr;
            for (Point* pta : ta->points) {
                if (pta != abc[0] && pta != abc[1] && pta != abc[2]) {
                    d = pta;
                    break;
                }
            }

            if (t->isPointInSphere(d)) {
                Vecteur v(*p, *d);
                bool intersect = intersect_droite_triangle(v, *p, *abc[0], *abc[1], *abc[2]);
                if (intersect) {
                    std::array<Tetraedre*, 3> tetras3 = flip23({t, ta});
                    ++flip23count;
                    tetras_to_check.push(tetras3[0]);
                    tetras_to_check.push(tetras3[1]);
                    tetras_to_check.push(tetras3[2]);
                    continue;
                }
                Tetraedre* tetra_adj_commun = tetra_adj_2_points(t, p, ta, d);
                if (!intersect && tetra_adj_commun != nullptr) {
                    flip32({t, ta, tetra_adj_commun});
                    ++flip32count;
                }
            }
        }
    }
}
