#include "triangulationdelaunay.h"

TriangulationDelaunay::TriangulationDelaunay(std::vector<Point>* tri_ret, std::vector<Point> points)
{

    float min_x = FLT_MAX, min_y = FLT_MAX, min_z = FLT_MAX;
    float max_x = FLT_MIN, max_y = FLT_MIN, max_z = FLT_MIN;
    for (Point point : points) {
        remaining_points.push(&point);
        if (point.x < min_x) min_x = point.x;
        else if (point.x > max_x) max_x = point.x;
        if (point.y < min_y) min_y = point.y;
        else if (point.y > max_y) max_y = point.y;
        if (point.z < min_z) min_z = point.z;
        else if (point.z > max_z) max_z = point.z;
    }

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


    Point p1(center_x, y_box_max, center_z);
    Point p2(x_box_min, y_box_min, z_box_min);
    Point p3(x_box_max, y_box_min, z_box_min);
    Point p4(center_x, y_box_min, z_box_max);

    Tetraedre* bigT = new Tetraedre(&p1, &p2, &p3, &p4);

    std::vector<Point> points_ret;
    tri_ret->push_back(p1);
    tri_ret->push_back(p2);
    tri_ret->push_back(p3);
    tri_ret->push_back(p4);

    std::vector<std::vector<float>> mat;
    std::vector<float> row1; row1.push_back(4); row1.push_back(2); row1.push_back(8); row1.push_back(4);
    std::vector<float> row2; row2.push_back(1); row2.push_back(4); row2.push_back(7); row2.push_back(9);
    std::vector<float> row3; row3.push_back(6); row3.push_back(5); row3.push_back(4); row3.push_back(7);
    std::vector<float> row4; row4.push_back(8); row4.push_back(8); row4.push_back(4); row4.push_back(3);
    mat.push_back(row1);mat.push_back(row2);mat.push_back(row3);mat.push_back(row4);

    tetraedres.push_front(bigT);

    triangulation();

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
        Vecteur v_test(*tetra->points[0], tetra->insphere_center);
        float dot1 = tetra->normales[0].dot_product(v_test);
        float dot2 = tetra->normales[1].dot_product(v_test);
        float dot3 = tetra->normales[2].dot_product(v_test);
        float dot4 = tetra->normales[3].dot_product(v_test);
        if (dot1 < 0 && dot2 < 0 && dot3 < 0 && dot4 < 0) {
            return tetra;
        }

        //on cherche la face vers laquelle se diriger ensuite
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
            //si les 3 points de la face appartiennent au triangle adjacent courant, c'est le bon
            for (int ind_points_tetra = 0 ; ind_points_tetra < 4 ; ++ind_points_tetra) {
                if (tetra_adj->points[ind_points_tetra] == p_face_1 ||
                        tetra_adj->points[ind_points_tetra] == p_face_2 ||
                        tetra_adj->points[ind_points_tetra] == p_face_3) {
                    ++common_points;
                }
            }
            if (common_points == 3) {
                tetra = tetra_adj;
            }
        }
    }


    return nullptr;
}

void TriangulationDelaunay::triangulation() {
    while (!remaining_points.empty()) {
        Point* cur_point = remaining_points.top();
        remaining_points.pop();
        Tetraedre* current_tetra = tetra_containing_point_walk(cur_point);
    }
}
