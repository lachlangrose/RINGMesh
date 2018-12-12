#include <ringmesh/ringmesh_tests_config.h>

#include <interpolation.h>

#include <geogram/mesh/mesh.h>
#include <geogram/basic/command_line.h>
//#include <ringmesh/io/geomodel/io_vtk.hpp>
#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>
#include <ringmesh/mesh/mesh_aabb.h>
#include <ringmesh/mesh/mesh_index.h>
#include <geogram/basic/attributes.h>

#include <ringmesh/mesh/volume_mesh.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ctime>
#include <iostream>
#include <vector>

using Eigen::MatrixXd;
using namespace RINGMesh;
Point3D::Point3D(double x, double y, double z) : _x(x), _y(y), _z(z) {}

InterfaceData::InterfaceData(double x, double y, double z, double v)
    : Point3D(x, y, z), _v(v) {}

DSI::DSI(const GeoModelMesh3D& mesh)
    : _mesh(mesh), _I(Eigen::MatrixXd::Zero(3, 4)), _c(0), _nc(0) {
    ;  // I(4,3);
    for (index_t c = 0; c < 3; c++) {
        _I(c, 0) = -1;
        _I(c, c + 1) = 1;
    }
}
void DSI::calculate_cell_gradient(index_t cell_global_index,
                                  Eigen::MatrixXd& grad,
                                  std::vector<index_t>& vert_idx) {
    Eigen::MatrixXd vertices(_mesh.cells.nb_vertices(cell_global_index), 3);

    for (index_t vi = 0; vi < _mesh.cells.nb_vertices(cell_global_index);
         vi++) {
        index_t vert =
            _mesh.cells.vertex(ElementLocalVertex(cell_global_index, vi));
        auto vertex = _mesh.cells.mesh().vertex(vert);
        for (index_t c = 0; c < 3; c++) {
            vertices(vi, c) = vertex[c];
        }
        vert_idx.push_back(vi);
    }
    Eigen::MatrixXd d_f(3, 3);
    for (index_t i = 0; i < 3; i++) {
        d_f(0, i) = vertices(i + 1, 0) - vertices(0, 0);
        d_f(1, i) = vertices(i + 1, 1) - vertices(0, 1);
        d_f(2, i) = vertices(i + 1, 2) - vertices(0, 2);
    }
    grad = d_f.inverse() * _I;
    return;
}
void DSI::add_constant_gradient(double w) {
    std::vector<double> A;
    std::vector<index_t> row;
    std::vector<index_t> col;
    std::vector<double> B;

    Eigen::MatrixXd d_n(3, 3);
    index_t Nc = 5;
    index_t Na = 4;
    int idc[Nc] = {0};  // std::vector<int> idc(Nc);        // constraint id
    double cstr[Nc] = {0.0};  // std::vector<double> cstr(Nc,0.0);

    int common_index;
    index_t next_available_position, position_to_write;
    for (index_t m = 0; m < _mesh.geomodel().nb_regions();
         m++) {  // loop over regions
        _nc = _mesh.cells.nb_vertices(m);
        for (index_t c = 0; c < _mesh.cells.nb_cells(m);
             c++) {  // loop over cells in region
            // get the cell gradient
            Eigen::MatrixXd e1(3, 3);
            std::vector<index_t> vert_idx;
            calculate_cell_gradient(c, e1, vert_idx);
            // get the neighbour gradient
            index_t cell_in_gmm = _mesh.cells.cell(m, c);
            for (index_t f = 0; f < _mesh.cells.nb_facets(cell_in_gmm);
                 f++) {  // loop over neighbours
                index_t n = _mesh.cells.adjacent(cell_in_gmm, f);
                if (n == GEO::NO_CELL) {
                    continue;
                }
                Eigen::MatrixXd shared_vertices(
                    _mesh.cells.nb_facet_vertices(
                        CellLocalFacet(cell_in_gmm, f)),
                    3);
                for (index_t f_v = 0; f_v < _mesh.cells.nb_facet_vertices(
                                                CellLocalFacet(cell_in_gmm, f));
                     f_v++) {  // loop over shared vertices
                    auto vert =
                        _mesh.cells.mesh().vertex(_mesh.cells.facet_vertex(
                            CellLocalFacet(cell_in_gmm, f), f_v));
                    for (index_t i = 0; i < 3; i++) {
                        shared_vertices(f_v, i) = vert[i];
                    }
                }  // loop
                Eigen::MatrixXd e2(3, 3);
                std::vector<index_t> n_vert_idx;
                calculate_cell_gradient(n, e2, n_vert_idx);
                Eigen::Vector3d v1;
                Eigen::Vector3d v2;
                // find the normal to the shared face
                for (index_t i = 0; i < 3; i++) {
                    v1(i) = shared_vertices(0, i) - shared_vertices(1, i);
                    v2(i) = shared_vertices(2, i) - shared_vertices(1, i);
                }
                auto norm = v1.cross(v2);
                // start with the original tetra
                for (index_t itr_left = 0; itr_left < Na; itr_left++) {
                    idc[itr_left] = vert_idx[itr_left];
                    for (index_t i = 0; i < 3; i++) {
                        cstr[itr_left] += norm(i) * e1(i, itr_left);
                    }
                }
                // process neighbour
                next_available_position = Na;
                for (index_t itr_right = 0; itr_right < Na; itr_right++) {

                    common_index = -1;
                    for (index_t itr_left = 0; itr_left < Na; itr_left++) {
                        if (idc[itr_left] == n_vert_idx[itr_right]) {
                            common_index = (int)itr_left;
                        }
                    }
                    position_to_write = 0;
                    if (common_index != -1) {
                        position_to_write = (index_t)common_index;
                    } else {
                        position_to_write = next_available_position;

                        next_available_position += 1;
                    }
                    idc[position_to_write] = n_vert_idx[itr_right];
                    for (index_t i = 0; i < 3; i++) {
                        cstr[position_to_write] -= norm(i) * e2(i, itr_right);
                    }
                }
                for (index_t i = 0; i < Nc; i++) {
                    // std::cout<<i<<" " <<cstr[i]<<" "<<idc[i]<<"
                    // "<<_c<<std::endl;
                    _A.push_back(cstr[i]);
                    _col.push_back(idc[i]);
                    _row.push_back(_c);
                    _B.push_back(0.0);
                    idc[i] = 0.0;
                    cstr[i] = 0;
                }
                // install constraint
            }      // end loop over neighbours
            _c++;  // next constraint adds a new row to the matrix
        }
    }
}

void DSI::add_control_points(std::vector<InterfaceData> interfacedata) {
    const VolumeAABBTree<3>& aabb3D = _mesh.cells.mesh().cell_aabb();
    Eigen::MatrixXd M(4, 4);
    for (const auto& iface : interfacedata) {

        M = Eigen::MatrixXd::Constant(4, 4, 1.0);
        vecn<3> point(iface.x(), iface.y(), iface.z());
        index_t containing_cell = aabb3D.containing_cell(point);
        if (containing_cell != GEO::NO_CELL) {
            for (index_t vi = 0; vi < _mesh.cells.nb_vertices(containing_cell);
                 vi++) {
                index_t vert_id =
                    _mesh.cells.vertex(ElementLocalVertex(containing_cell, vi));
                auto vert = _mesh.cells.mesh().vertex(vert_id);
                for (index_t j = 0; j < 3; j++) {
                    M(vi, j) = vert[j];
                }
            }
            Eigen::MatrixXd cp(1, 4);
            cp = Eigen::MatrixXd::Constant(1, 4, 1.0);
            cp(1) = iface.x();
            cp(2) = iface.y();
            cp(3) = iface.z();
            auto c = cp * M.inverse();
            for (index_t i = 0; i < _mesh.cells.nb_vertices(containing_cell);
                 i++) {
                _A.push_back(c(i));
                _row.push_back(_c);
                _col.push_back(
                    _mesh.cells.vertex(ElementLocalVertex(containing_cell, i)));
            }
            _B.push_back(iface._v);
            _c++;
        }
    }
}
void DSI::solve_system() {
    Eigen::SparseMatrix<double> A(_nc, _c);
    Eigen::VectorXd B(_nc);
    A.reserve(_A.size());
    for (index_t i = 0; i < _A.size(); i++) {
        A.insert(_col[i], _row[i]) = _A[i];
    }
    for (index_t i = 0; i < _nc; i++) {
        B(i) = _B[i];
    }
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > lscg;
    lscg.compute(A);

    auto x = lscg.solve(B);
    std::cout << "#iterations:     " << lscg.iterations() << std::endl;
    std::cout << "estimated error: " << lscg.error() << std::endl;
    for (index_t m = 0; m < _mesh.geomodel().nb_regions();
         m++) {  // loop over regions
        const Region3D& cur_reg = _mesh.geomodel().region(m);

        GEO::AttributesManager& reg_attr_mgr =
            cur_reg.vertex_attribute_manager();
        GEO::Attribute<double> interpolant_attr(reg_attr_mgr, "interpolant");

        for (index_t c = 0; c < _nc; c++) {
            interpolant_attr[c] = x(c);  //( rounded_volume % 2 == 0 );
        }
    }
}
void DSI::add_gradient_control_points(std::vector<PlanarData> planardata) {}
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
    std::vector<InterfaceData> ifaces;
    ifaces.push_back(InterfaceData(0.1, 0.1, 0.1, 0.0));
    ifaces.push_back(InterfaceData(5.1, 5.1, 5.1, 1.0));
    dsi->add_control_points(ifaces);
    dsi->solve_system();
    // double duration;
    // duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    //// Tetrahedralize the GeoModel
    //// tetrahedralize(geomodel, NO_ID, false);
    // std::string file_name_out(ringmesh_test_data_path);
     geomodel_save(geomodel, ringmesh_test_output_path + "geomodel" +
                                std::to_string(3) + "d." + "so");
    // index_t nb{0};
    // index_t nc{0};
    //// iterate over all regions in the geomodel
    // for (const auto& region : geomodel.regions()) {
    //    auto r = region.index();
    //    nb += region.nb_mesh_elements();
    //}
    // std::cout << "That took: " << duration << " for " << nb << " cells"
    //          << std::endl;

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
