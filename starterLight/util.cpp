#include "util.h"

int calcDeterminant(std::vector<std::vector<float>> Matrix) {
    //this function is written in c++ to calculate the determinant of matrix
    // it's a recursive function that can handle matrix of any dimension
    int det = 0; // the determinant value will be stored here
    if (Matrix.size() == 1)
    {
        return Matrix[0][0]; // no calculation needed
    }
    else if (Matrix.size() == 2)
    {
        //in this case we calculate the determinant of a 2-dimensional matrix in a
        //default procedure
        det = (Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0]);
        return det;
    }
    else
    {
        //in this case we calculate the determinant of a squared matrix that have
        // for example 3x3 order greater than 2
        for (int p = 0; p < Matrix[0].size(); p++)
        {
            //this loop iterate on each elements of the first row in the matrix.
            //at each element we cancel the row and column it exist in
            //and form a matrix from the rest of the elements in the matrix
            std::vector<std::vector<float>> TempMatrix; // to hold the shaped matrix;
            for (int i = 1; i < Matrix.size(); i++)
            {
                // iteration will start from row one cancelling the first row values
                std::vector<float> TempRow;
                for (int j = 0; j < Matrix[i].size(); j++)
                {
                    // iteration will pass all cells of the i row excluding the j
                    //value that match p column
                    if (j != p)
                    {
                        TempRow.push_back(Matrix[i][j]);//add current cell to TempRow
                    }
                }
                if (TempRow.size() > 0)
                    TempMatrix.push_back(TempRow);
                //after adding each row of the new matrix to the vector tempx
                //we add it to the vector temp which is the vector where the new
                //matrix will be formed
            }
            det = det + Matrix[0][p] * pow(-1, p) * calcDeterminant(TempMatrix);
            //then we calculate the value of determinant by using a recursive way
            //where we re-call the function by passing to it the new formed matrix
            //we keep doing this until we get our determinant
        }
        return det;
    }
}

bool intersect_droite_triangle(Vecteur& vec_direct_droite, Point& p_droite, Point& p1, Point& p2, Point& p3) {

    //on cherche point d'intersection entre droite et plan
    Vecteur normale_plan = Vecteur(p1,p2).cross_product(Vecteur(p1,p3));
    Point p_plan = p1;

    float d = -(normale_plan.dot_product(Vecteur(p_plan.x, p_plan.y, p_plan.y)));
    float haut = -(normale_plan.dot_product(Vecteur(p_droite.x, p_droite.y, p_droite.z)) + d);

    float bas = normale_plan.dot_product(vec_direct_droite);

    Vecteur v = vec_direct_droite * (haut / bas);
    Point p_inters = p_droite + Point(v.x, v.y, v.z);

    //on regarde si ce point d'intersection appartient au triangle
    Vecteur v_tri_1(p1,p2);
    Vecteur v_tri_2(p2,p3);
    Vecteur v_tri_3(p3,p1);

    Vecteur c1(p1, p_inters);
    Vecteur c2(p2, p_inters);
    Vecteur c3(p3, p_inters);

    if (normale_plan.dot_product(v_tri_1.cross_product(c1)) > 0 &&
        normale_plan.dot_product(v_tri_2.cross_product(c2)) > 0 &&
        normale_plan.dot_product(v_tri_3.cross_product(c3)) > 0) {
        return true;
    }
    return false;

}
