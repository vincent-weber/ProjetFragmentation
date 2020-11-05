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

    Tetraedre bigT(&p1, &p2, &p3, &p4);

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

    qDebug() << CalcDeterminant(mat);


    tetraedres.push_front(&bigT);

    triangulation();

}

Tetraedre* TriangulationDelaunay::tetra_containing_point_walk(Point* point) {
    //algorithme walk
    std::vector<Tetraedre*> visited_tetras;
    for (Tetraedre* tetra : tetraedres) {
        visited_tetras.push_back(tetra);

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
