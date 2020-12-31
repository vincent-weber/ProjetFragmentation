#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <math.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/internal/repair_extra.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh_CGAL;

typedef Mesh_CGAL::Vertex_index vertex_descriptor;
typedef Mesh_CGAL::Edge_index edge_descriptor;
typedef Mesh_CGAL::Face_index face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::Polygon_mesh_processing::parameters;

#include <stack>

#include "point.h"
#include "triangulationdelaunay.h"

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

struct BoiteEnglobante {
    float minX, maxX, minY, maxY, minZ, maxZ;

    BoiteEnglobante(float minX, float maxX, float minY, float maxY, float minZ, float maxZ) {
        this->minX = minX;
        this->maxX = maxX;
        this->minY = minY;
        this->maxY = maxY;
        this->minZ = minZ;
        this->maxZ = maxZ;
    }
};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    std::vector<Point*> listSeeds;
    int nbSeeds = 10;
    int gridSize = 5;
    int idPointImpact = 0;

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);
    BoiteEnglobante boiteEnglobante(MyMesh *_mesh);
    VertexHandle* find_vertex(Point& p, std::vector<VertexHandle>& handles);
    void deletesSeeds();

    Mesh_CGAL convert_open_mesh_to_cgal(MyMesh& openmesh_mesh);
    std::vector<Mesh_CGAL> convert_tetras_to_cgal(std::vector<Tetraedre>& tetras);
    void write_tetras_to_file(std::vector<Tetraedre>& tetras, std::string filename_pattern);
    void tests();

private slots:
    void on_pushButton_chargement_clicked();
    void on_pushButtonFragmentation_clicked();

    void on_spinBoxSeeds_valueChanged(int arg1);
    void on_spinBoxPointImpact_valueChanged(int arg1);
    void on_spinBoxGridSize_valueChanged(int arg1);

    void on_pushButtonRandom_clicked();
    void on_pushButtonImpact_clicked();
    void on_pushButtonGrid_clicked();

private:

    MyMesh mesh;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
