#include <ringmesh/ringmesh_tests_config.h>

#include <interpolation.h>

#include <geogram/basic/command_line.h>
//#include <ringmesh/io/geomodel/io_vtk.hpp>
#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/mesh_aabb.h>

#include <ringmesh/mesh/volume_mesh.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <ctime>

using Eigen::MatrixXd;
using namespace RINGMesh;
Point3D::Point3D(double x, double y, double z) : _x(x), _y(y), _z(z) {}

InterfaceData::InterfaceData(double x, double y, double z, double v)
    : Point3D(x, y, z), _v(v) {}

DSI::DSI(const GeoModelMesh3D& mesh) : _mesh(mesh) {}
void DSI::add_constant_gradient(double w) {
    Eigen::MatrixXd I(3, 4);

    Eigen::MatrixXd d_f(3, 3);
    Eigen::MatrixXd d_n(3, 3);
    I = Eigen::MatrixXd::Zero(3, 4);  // I(4,3);
    for (index_t c = 0; c < 3; c++) {
        I(c, 0) = -1;
        I(c, c + 1) = 1;
    }
    std::cout << I << std::endl;

    for (index_t m = 0; m < _mesh.geomodel().nb_regions(); m++) {
        for (index_t c = 0; c < _mesh.cells.nb_cells(m);
             c++) {  // It is also possible to iterate only on one element type
            index_t cell_in_gmm = _mesh.cells.cell(m, c);
            Eigen::MatrixXd vertices(_mesh.cells.nb_vertices(c), 3);
            for (index_t vi = 0; vi < _mesh.cells.nb_vertices(c); vi++) {
                index_t vert = _mesh.cells.vertex(ElementLocalVertex(c, vi));
                auto vertex = _mesh.cells.mesh().vertex(vert);
                for (index_t c = 0; c < 3; c++) {
                    vertices(vi, c) = vertex[c];
                }
                // d_f(0,0) =
                // std::cout<<"vert "<<vertex[0]<<" "<<vertex[1]<<"
                // "<<vertex[2]<<std::endl;
            }
            for (index_t i = 0; i < 3; i++) {
                d_f(0, i) = vertices(1, i) - vertices(0, i);
                d_f(1, i) = vertices(2, i) - vertices(0, i);
                d_f(2, i) = vertices(3, i) - vertices(0, i);
            }
            auto grad = d_f.inverse() * I;
            for (index_t f = 0; f < _mesh.cells.nb_facets(cell_in_gmm); f++) {
                index_t n = _mesh.cells.adjacent(
                    cell_in_gmm, f);  // Global cell adjacent index
            Eigen::MatrixXd vertices_n(_mesh.cells.nb_vertices(n), 3);
            for (index_t vi = 0; vi < _mesh.cells.nb_vertices(n); vi++) {
                index_t vert = _mesh.cells.vertex(ElementLocalVertex(n, vi));
                auto vertex = _mesh.cells.mesh().vertex(vert);
                for (index_t i = 0; i < 3; i++) {
                    vertices_n(vi, i) = vertex[i];
                }
                // d_f(0,0) =
                // std::cout<<"vert "<<vertex[0]<<" "<<vertex[1]<<"
                // "<<vertex[2]<<std::endl;
            }
            for (index_t c = 0; c < 3; c++) {
                d_f(0, i) = vertices_n(1, i) - vertices(0, 0);
                d_f(1, i) = vertices_n(2, i) - vertices(0, 1);
                d_f(2, i) = vertices_n(3, i) - vertices(0, 2);
            }
            }
        }
    }
}
void add_control_points(std::vector<InterfaceData> interfacedata) {}

void DSI::add_gradient_control_points(std::vector<PlanarData> planardata) {}
void DSI::calculate_constant_gradient(index_t c1, index_t c2,
                                      std::vector<int>& idc,
                                      std::vector<double>& c) {}
int main(int argc, char** argv) {
    using namespace RINGMesh;

    std::string file_name(ringmesh_test_data_path);
    file_name += argv[1];  //"modelA1.ml";

    // Check only model geometry
    GeoModel3D geomodel;
    bool loaded_model_is_valid = geomodel_load(geomodel, file_name);
    std::cout << "trying to make dsi" << std::endl;
    DSI* dsi = new DSI(geomodel.mesh);
    std::cout << "made a dsi" << std::endl;
    std::clock_t start;
    start = std::clock();
    dsi->add_constant_gradient();
    double duration;
    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    // Tetrahedralize the GeoModel
    // tetrahedralize(geomodel, NO_ID, false);
    std::string file_name_out(ringmesh_test_data_path);
    geomodel_save(geomodel, ringmesh_test_output_path + "geomodel" +
                                std::to_string(3) + "d." + "vtk");
    index_t nb{0};
    index_t nc{0};
    // iterate over all regions in the geomodel
    for (const auto& region : geomodel.regions()) {
        auto r = region.index();
        nb += region.nb_mesh_elements();
    }
    std::cout << "That took: " << duration << " for " << nb << " cells"
              << std::endl;

    //
    //        // Check validity of tetrahedralized model
    //        ValidityCheckMode checks{ ValidityCheckMode::GEOMETRY };
    //
    //        if( !is_geomodel_valid( geomodel, checks ) )
    //        {
    //            throw RINGMeshException( "RINGMesh Test",
    //                "Failed when tetrahedralize model ", geomodel.name(),
    //                ": the model becomes invalid." );
    //        }
}
