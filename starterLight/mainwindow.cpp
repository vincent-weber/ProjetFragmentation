#include "mainwindow.h"
#include "ui_mainwindow.h"

/* **** début de la partie boutons et IHM **** */

BoiteEnglobante MainWindow::boiteEnglobante(MyMesh *_mesh) {
    float minX = 0, maxX = 0, minY = 0, maxY = 0, minZ = 0, maxZ = 0;
    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++) {
        for (MyMesh::FaceVertexIter fv_it=mesh.fv_iter(*curFace); fv_it.is_valid(); ++fv_it) {
            if(_mesh->point(*fv_it)[0] < minX) minX = _mesh->point(*fv_it)[0];
            else if(_mesh->point(*fv_it)[0] > maxX) maxX = _mesh->point(*fv_it)[0];
            if(_mesh->point(*fv_it)[1] < minY) minY = _mesh->point(*fv_it)[1];
            else if(_mesh->point(*fv_it)[1] > maxY) maxY = _mesh->point(*fv_it)[1];
            if(_mesh->point(*fv_it)[2] < minZ) minZ = _mesh->point(*fv_it)[2];
            else if(_mesh->point(*fv_it)[2] > maxZ) maxZ = _mesh->point(*fv_it)[2];
        }
    }
    return BoiteEnglobante(minX, maxX, minY, maxY, minZ, maxZ);
}

float rand_float_between_two_values(float a, float b) {
    float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    return (b-a) * r + a;
}

std::vector<Point*> genererPointsDansBoite(MyMesh *_mesh, BoiteEnglobante& boite, int nb_points) {

    std::vector<Point*> points;
    for (int i = 0 ; i < nb_points ; ++i) {
        float x = rand_float_between_two_values(boite.minX, boite.maxX);
        float y = rand_float_between_two_values(boite.minY, boite.maxY);
        float z = rand_float_between_two_values(boite.minZ, boite.maxZ);

        Point* p_pile = new Point(x,y,z);
        points.push_back(p_pile);

        MyMesh::Point p;
        p[0] = x;
        p[1] = y;
        p[2] = z;
        VertexHandle vh = _mesh->add_vertex(p);
        _mesh->set_color(vh, MyMesh::Color(0, 0, 255));
        _mesh->data(vh).thickness = 5;
    }
    return points;
}

std::vector<Point*> genererPointsGrid(MyMesh *_mesh, BoiteEnglobante& boite, int size) {
    std::vector<Point*> points;
    float xSpace = (boite.maxX - boite.minX) / (size - 1);
    float ySpace = (boite.maxY - boite.minY) / (size - 1);
    float zSpace = (boite.maxZ - boite.minZ) / (size - 1);

    for (int x = 0 ; x < size ; x++) {
        for (int y = 0 ; y < size ; y++){
            for (int z = 0 ; z < size ; z++){
                float newX = boite.minX + xSpace * x + rand_float_between_two_values(-0.01, 0.01);
                float newY = boite.minY + ySpace * y + rand_float_between_two_values(-0.01, 0.01);
                float newZ = boite.minZ + zSpace * z + rand_float_between_two_values(-0.01, 0.01);
                Point* p_pile = new Point(newX,newY,newZ);
                points.push_back(p_pile);


                MyMesh::Point p;
                p[0] = newX;
                p[1] = newY;
                p[2] = newZ;
                VertexHandle vh = _mesh->add_vertex(p);
                _mesh->set_color(vh, MyMesh::Color(0, 0, 255));
                _mesh->data(vh).thickness = 5;
            }
        }
    }
    return points;
}

float distanceEuclidienne(MyMesh::Point pA, MyMesh::Point pB){
    return (float) sqrt(pow(pA[0]-pB[0], 2) + pow(pA[1]-pB[1], 2) + pow(pA[2]-pB[2], 2));
}

std::vector<Point*> genererPointsImpact(MyMesh *_mesh, BoiteEnglobante& boite, int nb_points, int idPoint) {
    std::vector<Point*> points;
    MyMesh::Point impact = _mesh->point(_mesh->vertex_handle(idPoint));
    float xDistMax = std::max(sqrt(pow(impact[0]-boite.maxX, 2)), sqrt(pow(impact[0]-boite.minX, 2)));
    float pasX = xDistMax / nb_points;
    float yDistMax = std::max(sqrt(pow(impact[1]-boite.maxY, 2)), sqrt(pow(impact[1]-boite.minY, 2)));
    float pasY = yDistMax / nb_points;
    float zDistMax = std::max(sqrt(pow(impact[2]-boite.maxZ, 2)), sqrt(pow(impact[2]-boite.minZ, 2)));
    float pasZ = zDistMax / nb_points;

    for (int i = 0 ; i < nb_points ; ++i) {
        float x = rand_float_between_two_values(impact[0]-(pasX*(i+1)), impact[0]+(pasX*(i+1)));
        float y = rand_float_between_two_values(impact[1]-(pasY*(i+1)), impact[1]+(pasY*(i+1)));
        float z = rand_float_between_two_values(impact[2]-(pasZ*(i+1)), impact[2]+(pasZ*(i+1)));

        if (x > boite.minX && y > boite.minY && z > boite.minZ && x < boite.maxX && y < boite.maxY && z < boite.maxZ){
            Point* p_pile = new Point(x,y,z);
            points.push_back(p_pile);

            MyMesh::Point p;
            p[0] = x;
            p[1] = y;
            p[2] = z;
            VertexHandle vh = _mesh->add_vertex(p);
            _mesh->set_color(vh, MyMesh::Color(0, 0, 255));
            _mesh->data(vh).thickness = 10;
        } else {
            i--;
        }
    }

    return points;
}

VertexHandle* MainWindow::find_vertex(Point& p, std::vector<VertexHandle>& handles) {
    for (VertexHandle& vh : handles) {
        MyMesh::Point ph = mesh.point(vh);
        if (p.x == ph[0] && p.y == ph[1] && p.z == ph[2]) {
            return &vh;
        }
    }
    return nullptr;
}

void MainWindow::write_tetras_to_file(std::vector<Tetraedre>& tetras, std::string filename_pattern) {
    unsigned file_number = 0;
    std::vector<MyMesh> meshes_tetra;
    std::vector<VertexHandle> handles;
    for (Tetraedre& tetra : tetras) {
        MyMesh mesh_tetra;
        for (int i = 0 ; i < 4 ; ++i) {
            Point* p = tetra.points[i];
            MyMesh::Point ph(p->x, p->y, p->z);
            VertexHandle vh = mesh_tetra.add_vertex(ph);
            mesh_tetra.set_color(vh, MyMesh::Color(0, 255, 0)); mesh.data(vh).thickness = 10;
            handles.push_back(vh);
        }
        for (int i = 0 ; i < 4 ; ++i) {
            bool is_ccw = tetra.orient_faces[i];
            if (!is_ccw) {
                mesh_tetra.add_face(handles[(i+1)%4], handles[i], handles[(i+2)%4]);
            }
            else {
                mesh_tetra.add_face(handles[i], handles[(i+1)%4], handles[(i+2)%4]);
            }
        }
        ++file_number;
        try
          {
            std::string filename("./" + filename_pattern + std::to_string(file_number) + ".obj");
            if ( !OpenMesh::IO::write_mesh(mesh_tetra, filename) )
            {
              qDebug() << "Cannot write mesh to file 'output.obj'";
              return;
            }
          }
          catch( std::exception& x )
          {
            qDebug() << x.what();
            return;
          }

        meshes_tetra.push_back(mesh_tetra);
        handles.clear();
    }
}

std::vector<Mesh_CGAL> MainWindow::convert_tetras_to_cgal(std::vector<Tetraedre>& tetras) {
    qDebug() << "Starting converting tetras.";
    std::vector<Mesh_CGAL> meshes_tetra;
    std::vector<vertex_descriptor> handles;
    for (Tetraedre& tetra : tetras) {
        Mesh_CGAL mesh_tetra;
        for (int i = 0 ; i < 4 ; ++i) {
            vertex_descriptor vd = mesh_tetra.add_vertex(Mesh_CGAL::Point(tetra.points[i]->x,tetra.points[i]->y,tetra.points[i]->z));
            handles.push_back(vd);
        }
        for (int i = 0 ; i < 4 ; ++i) {
            bool is_ccw = tetra.orient_faces[i];
            if (!is_ccw) {
                mesh_tetra.add_face(handles[(i+1)%4], handles[i], handles[(i+2)%4]);
            }
            else {
                mesh_tetra.add_face(handles[i], handles[(i+1)%4], handles[(i+2)%4]);
            }
        }
        meshes_tetra.push_back(mesh_tetra);
        handles.clear();
    }
    qDebug() << "Finished converting tetras.";
    return meshes_tetra;
}

// Finds if the given point p is already present in the mesh.
// If yes, returns the corresponding handle.
// Else, add the point to the mesh and return the corresponding handle.
vertex_descriptor find_or_add(Mesh_CGAL& mesh, Mesh_CGAL::Point p)
{
    Mesh_CGAL::Property_map<vertex_descriptor, K::Point_3> location = get(CGAL::vertex_point, mesh);
    Mesh_CGAL::Vertex_range::iterator  vb, ve;
    Mesh_CGAL::Vertex_range r = mesh.vertices();

    vb = boost::begin(r);
    vb = boost::end(r);

    for(boost::tie(vb, ve) = mesh.vertices(); vb != ve && location[*vb] != p; ++vb);
    if (vb != mesh.vertices_end())
        return *vb;
    else
        return mesh.add_vertex(p);
}

Mesh_CGAL MainWindow::convert_open_mesh_to_cgal(MyMesh& openmesh_mesh) {
    qDebug() << "Starting converting mesh to cgal";
    Mesh_CGAL final_mesh_cgal;
    for (MyMesh::FaceIter curFace = openmesh_mesh.faces_begin(); curFace != openmesh_mesh.faces_end(); curFace++) {
        std::vector<vertex_descriptor> vertices_face;
        vertices_face.reserve(3);
        for (MyMesh::FaceVertexIter fv_iter = openmesh_mesh.cfv_iter(*curFace); fv_iter.is_valid() ; ++fv_iter) {
            MyMesh::Point pom(openmesh_mesh.point(*fv_iter));
            Mesh_CGAL::Point p(pom[0], pom[1], pom[2]);
            vertices_face.emplace_back(find_or_add(final_mesh_cgal, p));
        }
        final_mesh_cgal.add_face(vertices_face[0], vertices_face[1], vertices_face[2]);
    }
    qDebug() << "Finished converting mesh to cgal";
    return final_mesh_cgal;
}

Mesh_CGAL compute_intersection(Mesh_CGAL mesh_1, Mesh_CGAL mesh_2, bool debug_info = false) {
    bool valid_intersection = PMP::corefine_and_compute_intersection(mesh_1, mesh_2, mesh_1, params::throw_on_self_intersection(false));
    bool removed = false;
    bool intersecting = false;

    if (valid_intersection) {
        intersecting = PMP::does_self_intersect(mesh_1, PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh_1)));
        if (intersecting) {
            removed = PMP::remove_self_intersections(mesh_1);
        }
    }

    if (debug_info) {
        if (!valid_intersection) {
            qDebug() << "Intersection failed.";
            //exit(0);
        } else {
            if (intersecting) {
                qDebug() << (intersecting ? "There are self-intersections in Intersection." : "There is no self-intersection in Intersection.");
                if (!removed) {
                    qDebug() << "Self intersections not removed.";
                } else {
                    qDebug() << "Self intersections removed.";
                }
            }
        }
    }
    return mesh_1;
}

bool save_mesh_cgal(Mesh_CGAL& mesh,std::string filename)
{
    // write mesh to output.obj
    try
    {
        std::ofstream output(filename);
        output << mesh;
    }
    catch( std::exception& x )
    {
        qDebug() << "Exception: " << x.what() << endl;
        return false;
    }
    return true;
}

// exemple pour charger un fichier .obj
void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    srand(time(NULL));

    //tests();exit(0);

    // on affiche le maillage
    displayMesh(&mesh);
}

void MainWindow::tests() {

    /*Point p1(1,1,1);
    Point p2(1,-1,-1);
    Point p3(-1,1,-1);
    Point p4(-1,-1,1);*/
    Point p1(-1,0,0);
    Point p2(1,0,0);
    Point p3(0,0,1);
    Point p4(0,1,0);
    Tetraedre tetra_test(&p1, &p2, &p3, &p4);
    std::vector<Point> points_tests;
    points_tests.push_back(Point(0,-0.5,0.2)); //TRUE
    points_tests.push_back(Point(0,-0.5,-0.2)); //FALSE
    points_tests.push_back(Point(0.1,-10,0.1)); //TRUE
    points_tests.push_back(Point(0.1,-10,-0.1)); //FALSE
    points_tests.push_back(Point(0,-0.1,10)); //FALSE
    /*for (int i = 0 ; i < 4 ; ++i) {
        qDebug() << "NORMALE " << i << " : " << tetra_test.normales[i].x << " - " << tetra_test.normales[i].y << " - " << tetra_test.normales[i].z;
    }
    qDebug() << "INSPHERE CENTER : " << tetra_test.insphere_center.x << " - " << tetra_test.insphere_center.y << " - " << tetra_test.insphere_center.z;*/
    /*NORMALE  0  :  0  -  -1  -  0
    NORMALE  1  :  1  -  0  -  0
    NORMALE  2  :  -1  -  1  -  1
    NORMALE  3  :  0  -  0  -  -1*/

    for (Point point_test : points_tests) {
        Vecteur vdir(point_test, p4);
        qDebug() << intersect_droite_triangle(vdir, p4, p1,p2,p3);
        /*Vecteur v_test(*tetra_test.points[0], point_test);
        Vecteur v_test2(*tetra_test.points[1], point_test);
        float dot1 = tetra_test.normales[0].dot_product(v_test);
        float dot2 = tetra_test.normales[1].dot_product(v_test2);
        float dot3 = tetra_test.normales[2].dot_product(v_test);
        float dot4 = tetra_test.normales[3].dot_product(v_test);
        if (dot1 < 0 && dot2 < 0 && dot3 < 0 && dot4 < 0) {
            //qDebug() << "POINT DANS TETRAEDRE";
        }
        else {
            //qDebug() << "POINT HORS TETRAEDRE";
        }
        if (tetra_test.isPointInSphere(&point_test)) {
            //qDebug() << "POINT DANS SPHERE";
        } else {
            //qDebug() << "POINT HORS SPHERE";
        }*/
    }

    qDebug() << "CENTRE CERCLE CIRONSCRIT LOL : " << tetra_test.circumsphere_center;


    /*MyMesh::Point p_1; p_1[0] = p1.x; p_1[1] = p1.y; p_1[2] = p1.z; VertexHandle vh1 = mesh.add_vertex(p_1);
    mesh.set_color(vh1, MyMesh::Color(0, 255, 0)); mesh.data(vh1).thickness = 10;
    MyMesh::Point p_2; p_2[0] = p2.x; p_2[1] = p2.y; p_2[2] = p2.z; VertexHandle vh2 = mesh.add_vertex(p_2);
    mesh.set_color(vh2, MyMesh::Color(0, 255, 0)); mesh.data(vh2).thickness = 10;
    MyMesh::Point p_3; p_3[0] = p3.x; p_3[1] = p3.y; p_3[2] = p3.z; VertexHandle vh3 = mesh.add_vertex(p_3);
    mesh.set_color(vh3, MyMesh::Color(0, 255, 0)); mesh.data(vh3).thickness = 10;
    MyMesh::Point p_4; p_4[0] = p4.x; p_4[1] = p4.y; p_4[2] = p4.z; VertexHandle vh4 = mesh.add_vertex(p_4);
    mesh.set_color(vh4, MyMesh::Color(0, 255, 0)); mesh.data(vh4).thickness = 10;

    MyMesh::Point ps; ps[0] = tetra_test.insphere_center.x; ps[1] = tetra_test.insphere_center.y; ps[2] = tetra_test.insphere_center.z;
    VertexHandle vhs = mesh.add_vertex(ps);
    mesh.set_color(vhs, MyMesh::Color(0, 0, 255)); mesh.data(vhs).thickness = 10;

    displayMesh(&mesh);
    return;*/

    /*Point pb(1, 1, 1);
    Point pc(1,-1,-1);
    Point pd(-1,1,-1);
    Point pe(-1,-1,1);
    Tetraedre tetra2(&pc, &pb, &pe, &pd);
    Point ptest(3,0,0);
    qDebug() << "TEST IS POINT IN SPHERE";
    qDebug() << tetra2.isPointInSphere(&ptest);
    return;*/

}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButtonFragmentation_clicked()
{
    TriangulationDelaunay td(listSeeds, false);
    //return;
    std::vector<Tetraedre> tetras_td;
    for (Tetraedre* tetra : td.tetraedres) {
        tetras_td.push_back(*tetra);
    }

    Mesh_CGAL mesh_cgal = convert_open_mesh_to_cgal(mesh);
    std::vector<Mesh_CGAL> tetras_cgal = convert_tetras_to_cgal(tetras_td);

    std::string filename = "fragment";
    unsigned cnt = 1;
    qDebug() << "Starting computing fragments";
    for (Mesh_CGAL& tetra_cgal : tetras_cgal) {
        Mesh_CGAL mesh_intersect = compute_intersection(tetra_cgal, mesh_cgal, true);
        save_mesh_cgal(mesh_intersect, "./" + filename + std::to_string(cnt++) + ".off");
    }
    qDebug() << "Finished computing fragments";

    return;

    write_tetras_to_file(tetras_td, "tetra");
    write_tetras_to_file(td.tetras_debug, "tetra_debug");
}

void MainWindow::on_spinBoxSeeds_valueChanged(int arg1)
{
    nbSeeds = arg1;
}


void MainWindow::on_spinBoxPointImpact_valueChanged(int arg1)
{
    if(arg1 < (int)mesh.n_vertices()){
        idPointImpact = arg1;
        resetAllColorsAndThickness(&mesh);
        mesh.set_color(mesh.vertex_handle(idPointImpact), MyMesh::Color(255,0,0));
        mesh.data(mesh.vertex_handle(idPointImpact)).thickness = 5;
        displayMesh(&mesh);
    }
}

void MainWindow::on_spinBoxGridSize_valueChanged(int arg1)
{
    gridSize = arg1;
    displayMesh(&mesh);
}

void MainWindow::deletesSeeds(){
    for(int i = 1; i <= listSeeds.size(); i++){
        mesh.delete_vertex(mesh.vertex_handle(mesh.n_vertices()-i), false);
    }
    mesh.garbage_collection();
    listSeeds.clear();
}

void MainWindow::on_pushButtonRandom_clicked()
{
    deletesSeeds();
    BoiteEnglobante boite_englobante = boiteEnglobante(&mesh);
    listSeeds = genererPointsDansBoite(&mesh, boite_englobante, nbSeeds);
    displayMesh(&mesh);
}

void MainWindow::on_pushButtonImpact_clicked()
{
    deletesSeeds();
    BoiteEnglobante boite_englobante = boiteEnglobante(&mesh);
    listSeeds = genererPointsImpact(&mesh, boite_englobante, nbSeeds, idPointImpact);
    displayMesh(&mesh);
}

void MainWindow::on_pushButtonGrid_clicked()
{
    deletesSeeds();
    BoiteEnglobante boite_englobante = boiteEnglobante(&mesh);
    listSeeds = genererPointsGrid(&mesh, boite_englobante, gridSize);
    displayMesh(&mesh);
}
