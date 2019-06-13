#include <interpolation.h>
#include <mesher.h>
#include <pybind11/pybind11.h>
#include <ringmesh/ringmesh_tests_config.h>

#include <geogram/basic/attributes.h>
#include <geogram/basic/command_line.h>
#include <geogram/mesh/mesh.h>
#include <chrono>
//#include <ringmesh/io/geomodel/io_vtk.hpp>
#include <geogram/basic/attributes.h>
#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>
#include <ringmesh/mesh/mesh_aabb.h>
#include <ringmesh/mesh/mesh_index.h>

#include <ringmesh/mesh/volume_mesh.h>

#include <math.h>
#include <Eigen/SPQRSupport>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iterator>

#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#define DEG2RAD 0.0174533

using Eigen::MatrixXd;
using namespace RINGMesh;
using namespace RINGMesh;

BoundingBoxMesher::BoundingBoxMesher(GeoModel3D& geomodel)
    : geomodel_(geomodel) {}

void BoundingBoxMesher::mesh(std::vector<std::vector<double> > bounding_box,
			     double max_volume) {
	GeoModelBuilder3D builder(geomodel_);
	vec3 v0(bounding_box[0][0], bounding_box[1][0], bounding_box[2][0]);
	vec3 v1(bounding_box[0][1], bounding_box[1][0], bounding_box[2][0]);
	vec3 v2(bounding_box[0][1], bounding_box[1][1], bounding_box[2][0]);
	vec3 v3(bounding_box[0][0], bounding_box[1][1], bounding_box[2][0]);
	vec3 v4(bounding_box[0][0], bounding_box[1][0], bounding_box[2][1]);
	vec3 v5(bounding_box[0][1], bounding_box[1][0], bounding_box[2][1]);
	vec3 v6(bounding_box[0][1], bounding_box[1][1], bounding_box[2][1]);
	vec3 v7(bounding_box[0][0], bounding_box[1][1], bounding_box[2][1]);
	std::vector<index_t> triangles;
	triangles.push_back(0);
	triangles.push_back(1);
	triangles.push_back(2);
	triangles.push_back(0);
	triangles.push_back(2);
	triangles.push_back(3);

	std::vector<index_t> surface_facet_ptr;
	surface_facet_ptr.push_back(0);
	surface_facet_ptr.push_back(3);
	surface_facet_ptr.push_back(6);

	std::vector<vec3> vertices(4);
	gmme_id id;

	id = builder.topology.create_mesh_entity(Surface3D::type_name_static());
	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v2;
	vertices[3] = v3;
	builder.geometry.set_surface_geometry(id.index(), vertices, triangles,
					      surface_facet_ptr);

	id = builder.topology.create_mesh_entity(Surface3D::type_name_static());
	vertices[0] = v1;
	vertices[1] = v5;
	vertices[2] = v6;
	vertices[3] = v2;
	builder.geometry.set_surface_geometry(id.index(), vertices, triangles,
					      surface_facet_ptr);

	id = builder.topology.create_mesh_entity(Surface3D::type_name_static());
	vertices[0] = v4;
	vertices[1] = v5;
	vertices[2] = v6;
	vertices[3] = v7;
	builder.geometry.set_surface_geometry(id.index(), vertices, triangles,
					      surface_facet_ptr);

	id = builder.topology.create_mesh_entity(Surface3D::type_name_static());
	vertices[0] = v0;
	vertices[1] = v3;
	vertices[2] = v7;
	vertices[3] = v4;
	builder.geometry.set_surface_geometry(id.index(), vertices, triangles,
					      surface_facet_ptr);

	id = builder.topology.create_mesh_entity(Surface3D::type_name_static());
	vertices[0] = v3;
	vertices[1] = v2;
	vertices[2] = v6;
	vertices[3] = v7;
	builder.geometry.set_surface_geometry(id.index(), vertices, triangles,
					      surface_facet_ptr);

	id = builder.topology.create_mesh_entity(Surface3D::type_name_static());
	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v5;
	vertices[3] = v4;
	builder.geometry.set_surface_geometry(id.index(), vertices, triangles,
					      surface_facet_ptr);
	builder.build_lines_and_corners_from_surfaces();
	builder.build_regions_from_lines_and_surfaces();
	builder.end_geomodel();

	//// Tetrahedralize the GeoModel
	tetrahedralize(geomodel_, NO_ID, true, max_volume);
}

bool CellNeighbourIndex::operator==(const CellNeighbourIndex& rhs) const {
	if (rhs._cell1 == _cell1) {
		return rhs._cell2 == _cell2;
	}
	if (rhs._cell2 == _cell1) {
		return rhs._cell1 == _cell2;
	}
	return false;
}

bool CellNeighbourIndex::operator!=(const CellNeighbourIndex& rhs) const {
	return !operator==(rhs);
}
Point3D::Point3D(double x, double y, double z) : _x(x), _y(y), _z(z) {}

InterfaceData::InterfaceData(double x, double y, double z, double v)
    : Point3D(x, y, z), _v(v) {}
PlanarData::PlanarData(double x, double y, double z, Eigen::Vector3d v)
    : Point3D(x, y, z), _v(v / v.norm()) {}
PlanarData::PlanarData(double x, double y, double z, double strike, double dip)
    : Point3D(x, y, z) {
	double vx =
	    cos(-1.0 * (strike * DEG2RAD)) * cos(-1.0 * (dip * DEG2RAD));
	double vy =
	    sin(-1.0 * (strike * DEG2RAD)) * cos(-1.0 * (dip * DEG2RAD));
	double vz = sin(-1.0 * (dip * DEG2RAD));
	// Get strike vector - vp
	double vpx = -1.0 * vy;
	double vpy = vx;
	double vpz = 0;
	double Nx = (vy * vpz) - (vz * vpy);
	double Ny = (vz * vpx) - (vx * vpz);
	double Nz = (vx * vpy) - (vy * vpx);
	_v = Eigen::Vector3d(Nx, Ny, Nz);
	_v /= _v.norm();
}
double PlanarData::operator[](int x) {}
DSI::DSI(const GeoModelMesh3D& mesh)
    : _mesh(mesh), _I(Eigen::MatrixXd::Zero(3, 4)), _c(0), _nc(0) {
	;  // I(4,3);
	for (index_t c = 0; c < 3; c++) {
		_I(c, 0) = -1;
		_I(c, c + 1) = 1;
	}
	Logger::out("DSI", "Using DSI interpolation on Volumetric Mesh");
}
void DSI::calculate_cell_gradient(index_t cell_global_index,
				  Eigen::MatrixXd& grad,
				  std::vector<index_t>& vert_idx) {
	Eigen::MatrixXd vertices(
	    int(_mesh.cells.nb_vertices(cell_global_index)), 3);
	// vertices
	// x0 y0 z0 //
	// x1 y1 z1 //
	// x2 y2 z2 //
	// x3 y3 z3 //
	// x4 y4 z4 //
	for (index_t vi = 0; vi < _mesh.cells.nb_vertices(cell_global_index);
	     vi++) {
		index_t vert = _mesh.cells.vertex(
		    ElementLocalVertex(cell_global_index, vi));
		auto vertex = _mesh.cells.mesh().vertex(vert);
		for (index_t c = 0; c < 3; c++) {
			vertices(vi, c) = vertex[c];
		}
		vert_idx.push_back(vert);
	}
	Eigen::MatrixXd m(3, 3);  // from A.4 frank
	// x1-x0 y1-y0 z1-z0 //
	// x2-x0 y2-y0 z2-z0 //
	// x3-x0 y3-y0 z3-z0 //

	for (index_t i = 0; i < 3; i++) {
		m(i, 0) = vertices(i + 1, 0) - vertices(0, 0);
		m(i, 1) = vertices(i + 1, 1) - vertices(0, 1);
		m(i, 2) = vertices(i + 1, 2) - vertices(0, 2);
	}
	grad = m.inverse() * _I;
	return;
}
void DSI::add_constant_gradient(double w) {
	Logger::out("DSI", "Adding constant gradient");
	std::vector<double> A;
	std::vector<index_t> row;
	std::vector<index_t> col;
	std::vector<double> B;

	Eigen::MatrixXd d_n(3, 3);
	index_t Nc = 5;
	index_t Na = 4;
	int idc[Nc] = {0};  // std::vector<int> idc(Nc);        // constraint id
	double cstr[Nc] = {0.0};  // std::vector<double> cstr(Nc,0.0);
	std::vector<CellNeighbourIndex> visited;
	int common_index;
	index_t next_available_position, position_to_write;
	bool isVisited = false;
	for (index_t m = 0; m < _mesh.geomodel().nb_regions();
	     m++) {  // loop over regions
		_nc = _mesh.geomodel().region(m).nb_vertices();
		for (index_t c = 0; c < _mesh.cells.nb_cells(m);
		     c++) {  // loop over cells in region
			// get the cell gradient
			Eigen::MatrixXd e1(3, 3);
			std::vector<index_t> vert_idx;
			calculate_cell_gradient(c, e1, vert_idx);
			// get the index of the cell to finds its neighbour for each face
			index_t cell_in_gmm = _mesh.cells.cell(m, c);
			for (index_t f = 0;
			     f < _mesh.cells.nb_facets(cell_in_gmm);
			     f++) {  // loop over neighbours
				index_t n =
				    _mesh.cells.adjacent(cell_in_gmm, f);
				//check whether this pair of cells has already been added to the matrix
				CellNeighbourIndex pair(c, n);
				isVisited = false;
				for (const auto& it : visited) {
					if (it == pair) {
						isVisited = true;
						break;
					}
				}
				if (isVisited == true) {
					continue;
				}
				//if not add it to the container
				visited.push_back(pair);
				if (n == GEO::NO_CELL) {
					continue;
				}  // endif no cell
				Eigen::MatrixXd shared_vertices(
				    int(_mesh.cells.nb_facet_vertices(
					CellLocalFacet(cell_in_gmm, f))),
				    3);
				for (index_t f_v = 0;
				     f_v < _mesh.cells.nb_facet_vertices(
					       CellLocalFacet(cell_in_gmm, f));
				     f_v++) {  // loop over shared vertices
					auto vert = _mesh.cells.mesh().vertex(
					    _mesh.cells.facet_vertex(
						CellLocalFacet(cell_in_gmm, f),
						f_v));
					for (index_t i = 0; i < 3; i++) {
						shared_vertices(f_v, i) =
						    vert[i];
					}
				}  // loop faces
				Eigen::MatrixXd e2(3, 3);
				std::vector<index_t> n_vert_idx;
				calculate_cell_gradient(n, e2, n_vert_idx);
				Eigen::Vector3d v1;
				Eigen::Vector3d v2;
				// find the normal to the shared face
				for (index_t i = 0; i < 3; i++) {
					v1(i) = shared_vertices(0, i) -
						shared_vertices(1, i);
					v2(i) = shared_vertices(2, i) -
						shared_vertices(1, i);
				}
				auto norm = v1.cross(v2);
				// start with the original tetra
				// process neighbour
				//next_available_position = Na;
				for (index_t itr_right = 0; itr_right < Na;
				     itr_right++) {
					common_index = -1;
					for (index_t itr_left = 0;
					     itr_left < Na; itr_left++) {
						if (idc[itr_left] ==
						    n_vert_idx[itr_right]) {
							common_index =
							    (int)itr_left;
						}
					}
					position_to_write = 0;
					if (common_index != -1) {
						position_to_write =
						    (index_t)common_index;
					} else {
						position_to_write =
						    next_available_position;

						next_available_position += 1;
					}
					idc[position_to_write] =
					    n_vert_idx[itr_right];
					for (index_t i = 0; i < 3; i++) {
						cstr[position_to_write] -=
						    norm(i) * e2(i, itr_right);
					}
				}
				for (index_t i = 0; i < Nc; i++) {
					_A.push_back(cstr[i] * w);
					_col.push_back(idc[i]);
					_row.push_back(_c);
					idc[i] = 0.0;
					cstr[i] = 0;
				}
				_c++;  // next constraint adds a new row to the
				       // matrix
				_B.push_back(0.0);
				// install constraint
			}  // end loop over neighbours
		}
	}
}

void DSI::add_control_points(std::vector<InterfaceData> interfacedata,
			     double w) {
	const VolumeAABBTree<3>& aabb3D = _mesh.cells.mesh().cell_aabb();
	Eigen::MatrixXd M(4, 4);
	for (const auto& iface : interfacedata) {
		M = MatrixXd::Constant(4, 4, 1.0);
		vecn<3> point(iface.x(), iface.y(), iface.z());
		index_t containing_cell = aabb3D.containing_cell(point);
		if (containing_cell != GEO::NO_CELL) {
			auto a = _mesh.cells.mesh().vertex(
			     _mesh.cells.vertex(
				ElementLocalVertex(containing_cell, 0)));
			auto b = _mesh.cells.mesh().vertex(
			     _mesh.cells.vertex(
				ElementLocalVertex(containing_cell, 1)));
			auto c = _mesh.cells.mesh().vertex(
			     _mesh.cells.vertex(
				ElementLocalVertex(containing_cell, 2)));
			auto d = _mesh.cells.mesh().vertex(
			     _mesh.cells.vertex(
				ElementLocalVertex(containing_cell, 3)));

			auto vap = point - a;
			auto vbp = point - b;
			auto vcp = point - c;
			auto vdp = point - d;
			auto vab = b - a;
			auto vac = c - a;
			auto vad = d - a;
			auto vbc = c - b;
			auto vbd = d - b;
			
			vecn<4> bc;
			bc[0] =  GEO::dot(vbp,GEO::cross(vbd, vbc)) / 6.;
			bc[1] =  GEO::dot(vap,GEO::cross(vac, vad)) / 6.;
			bc[2] =  GEO::dot(vap,GEO::cross(vad, vab)) / 6.;
			bc[3] =  GEO::dot(vap,GEO::cross(vab, vac)) / 6.;
			double v = GEO::dot(vab,GEO::cross(vac,vad)) / 6.;
			bc/=v;

			//vecn<3> vb =  (vbd ^ vbc);
			//vecn<3> vc =  (vbd ^ vbc);
			//vecn<3> vd =  (vbd ^ vbc);
			//for (index_t vi = 0;
			//     vi < _mesh.cells.nb_vertices(containing_cell);
			//     vi++) {
			//	for (index_t j = 0; j < 3; j++) {
			//		M(vi, j) = vert[j];
			//	}
			//}
			//Eigen::MatrixXd cp(1, 4);
			//cp = Eigen::MatrixXd::Constant(1, 4, 1.0);
			//cp(1) = iface.x();
			//cp(2) = iface.y();
			//cp(3) = iface.z();
			//auto c = cp * M.inverse();
			for (index_t i = 0;
			     i < _mesh.cells.nb_vertices(containing_cell);
			     i++) {
				_A.push_back(bc[i]);
				_row.push_back(_c);
				_col.push_back(_mesh.cells.vertex(
				    ElementLocalVertex(containing_cell, i)));
			}
			_B.push_back(iface._v);
			_c++;
		}
	}
}

void DSI::add_gradient_control_points(std::vector<PlanarData> planardata,
				      double w) {
	const VolumeAABBTree<3>& aabb3D = _mesh.cells.mesh().cell_aabb();

	Eigen::MatrixXd e(3, 3);
	Eigen::VectorXd scalar_product(4);
	for (const auto& planar : planardata) {
		double length;
		vecn<3> point(planar.x(), planar.y(), planar.z());
		Logger::out("DSI",
			    "Adding gradient control point: ", planar._v);
		index_t containing_cell = aabb3D.containing_cell(point);
		if (containing_cell != GEO::NO_CELL) {
			std::vector<index_t> vert_idx;
			calculate_cell_gradient(containing_cell, e, vert_idx);
			for (index_t j = 0; j < 3; j++) {
				length = 0;
				for (index_t i = 0; i < 4; i++) {  // n nodes
					length += e(j, i) * e(j, i);
				}

				length = sqrt(length);
				for (index_t i = 0; i < _mesh.cells.nb_vertices(
							    containing_cell);
				     i++) {
					_A.push_back(e(j, i) / length);
					_row.push_back(_c);
					_col.push_back(_mesh.cells.vertex(
					    ElementLocalVertex(containing_cell,
							       i)));
				}
				double g;
				if (j == 0) g = planar._v.x();
				if (j == 1) g = planar._v.y();
				if (j == 2) g = planar._v.z();
				_B.push_back(g);
				_c++;
			}
		}
	}
}
void DSI::solve_system(int solver) {
	int n = 4;
	Eigen::setNbThreads(n);
	Eigen::SparseMatrix<double> A(_c, _nc);
	Eigen::VectorXd B(_c);
	A.reserve(_A.size());
	for (index_t i = 0; i < _A.size(); i++) {
		A.insert(_row[i], _col[i]) = _A[i];
	}
	Logger::out("DSI", "Interpolation matrix has ", A.cols(),
		    " columns and ", A.rows(), " rows");
	for (index_t i = 0; i < _c; i++) {
		B(i) = _B[i];
	}
	Eigen::VectorXd x;
	switch (solver) {
		case 0: {
			Logger::out("DSI",
				    "Using Eigen Least Squares Conjugate "
				    "Gradient solver");
			Eigen::LeastSquaresConjugateGradient<
			    Eigen::SparseMatrix<double> >
			    lscg;
			// lscg.setTolerance(0.001);
			lscg.compute(A);
			if (lscg.info() != Eigen::Success) {
				Logger::err("DSI",
					    "LSQR solver failed to initialise");
			}

			x = lscg.solve(B);
			if (lscg.info() != Eigen::Success) {
				if (lscg.info() == Eigen::NoConvergence) {
					Logger::err(
					    "DSI",
					    "LSQR solver did not converge");
				} else {
					Logger::err(
					    "DSI",
					    "LSQR solver failed to solve");
				}
			}
			Logger::out("DSI", "Solver took: ", lscg.iterations(),
				    "iterations ");
			Logger::out(" DSI",
				    " Estimated solver error: ", lscg.error());
			;
			break;
		}
		case 1: {
			Logger::out("DSI", "Using SuiteSparse SPQR solver");
			Eigen::SPQR<Eigen::SparseMatrix<double> > spqrsolver(A);
			x = spqrsolver.solve(B);
			break;
		}
		default: {
			Logger::err("DSI",
				    "Undefined flag for solver try 0 for Eigen "
				    "or 1 for SPQR");
			return;
		}
	}

	for (index_t m = 0; m < _mesh.geomodel().nb_regions();
	     m++) {  // loop over regions
		const Region3D& cur_reg = _mesh.geomodel().region(m);

		GEO::AttributesManager& reg_attr_mgr =
		    cur_reg.vertex_attribute_manager();
		GEO::Attribute<double> interpolant_attr(reg_attr_mgr,
							"interpolant");
		for (index_t c = 0; c < _nc; c++) {
			interpolant_attr[c] =
			    x(c);  //( rounded_volume % 2 == 0 );
		}
	}
}
int main(int argc, char** argv) {
	using namespace RINGMesh;

	std::string file_name(ringmesh_test_data_path);
	file_name += argv[1];
	// file_name += argv[1];  //"modelA1.ml";
	std::ifstream in_file(file_name);
	if (in_file.good() == false) {
		return 0;
	}
	std::string string;
	int flag = 0;
	std::vector<std::vector<double> > bb;

	std::vector<PlanarData> pdata;
	std::vector<InterfaceData> ifaces;
	while (in_file.good()) {
		std::getline(in_file, string);
		if (string.empty()) {
			continue;
		}
		if (string.compare("BB") == 0) {
			flag = 1;
			continue;
		}
		if (string.compare("GCP") == 0) {
			flag = 2;
			continue;
		}
		if (string.compare("CP") == 0) {
			flag = 3;
			continue;
		}
		std::istringstream buffer(string);
		switch (flag) {
			case 1: {
				std::vector<double> line(
				    (std::istream_iterator<double>(buffer)),
				    std::istream_iterator<double>());
				std::vector<double> tmp;
				for (std::vector<double>::iterator it =
					 line.begin();
				     it != line.end(); it++) {
					tmp.push_back(*it);
				}
				bb.push_back(tmp);
				continue;
			}
			case 2: {
				std::vector<double> line(
				    (std::istream_iterator<double>(buffer)),
				    std::istream_iterator<double>());
				if (line.size() == 5) {
					pdata.push_back(PlanarData(
					    line[0], line[1], line[2], line[3],
					    line[4]));
				}
				if (line.size() == 6) {
					pdata.push_back(PlanarData(
					    line[0], line[1], line[2],
					    Eigen::Vector3d(line[3], line[4],
							    line[5])));
				}

				continue;
			}
			case 3: {
				std::vector<double> line(
				    (std::istream_iterator<double>(buffer)),
				    std::istream_iterator<double>());
				ifaces.push_back(InterfaceData(
				    line[0], line[1], line[2], line[3]));

				continue;
			}
			default:
				continue;
		}
	}

	// Check only model geometry
	// GeoModel3D geomodel;
	GeoModel3D geomodel;
	// bb.push_back({0., 10.});
	// bb.push_back({0., 10.});
	// bb.push_back({-5., 1.});
	BoundingBoxMesher mesher(geomodel);
	mesher.mesh(bb, atof(argv[2]));
	// bool loaded_model_is_valid = geomodel_load(geomodel,
	// file_name);
	DSI* dsi = new DSI(geomodel.mesh);
	auto start = std::chrono::high_resolution_clock::now();
	dsi->add_constant_gradient();
	auto duration = (std::chrono::high_resolution_clock::now() - start);
	Logger::out(
	    "DSI", "Assembling CG constraint took: ",
	    std::chrono::duration_cast<std::chrono::seconds>(duration).count(),
	    " seconds");
	// ifaces.push_back(InterfaceData(1.1, 1.1, 1.1, 0.0));

	// ifaces.push_back(InterfaceData(2.5,3.52,0.4, 0.0));
	// ifaces.push_back(InterfaceData(8.1,9.1,3.1, 0.0));
	// ifaces.push_back(InterfaceData(7.5,1.5,0.4, 0.0));
	// fault_surface.add_point([1.1,1.1,0.1],0)
	//#fault_surface.add_point([8.1,9.1,3.1],10)
	// fault_surface.add_point([2.5,3.52,0.4],0)
	// fault_surface.add_point([5.3,2.2,0.4],0)
	// fault_surface.add_point([7.5,1.5,0.4],0)
	//
	//
	// fault_surface.add_strike_and_dip([2.5,3.52,0.4],0.,0.)
	// fault_surface.add_strike_and_dip([5.,2.52,0.4],45.,50.)
	// fault_surface.add_strike_and_dip([5.3,2.2,0.4],45.,50.)
	// ifaces.push_back(InterfaceData(0, 0, 0, 0.0));

	// ifaces.push_back(InterfaceData(2.5, 0.5, 0, 0.0));
	// ifaces.push_back(InterfaceData(5, 2.5, 0, 0.0));
	// ifaces.push_back(InterfaceData(9, 4, 0, 0.0));

	// pdata.push_back(PlanarData(2.5, .5, 0, 90.,90.));
	////pdata.push_back(PlanarData(5, 2.5, 0, 135.,90.));
	// pdata.push_back(PlanarData(9, 4, 0, 90.,90.));
	// pdata.push_back(PlanarData(9, 4, 0, 190.,80.));
	// pdata.push_back(PlanarData(9, 4, 0, 190.,80.));

	// pdata.push_back(PlanarData(-5, -5, 0, Eigen::Vector3d(.1,
	// -.5, .8))); for (double x = -10.0; x < 10.0; x += 1.0) {
	//	ifaces.push_back(InterfaceData(x, 0.1, -9.1, 0.0));
	//	ifaces.push_back(InterfaceData(x, 0.1, 0., 0.0));

	//	ifaces.push_back(InterfaceData(x, 5.1, -5.1, 9.0));
	//	ifaces.push_back(InterfaceData(x, 5.1, 0., 9.0));
	//}
	dsi->add_control_points(ifaces);
	dsi->add_gradient_control_points(pdata, 10);
	auto temp = (std::chrono::high_resolution_clock::now() - start);

	Logger::out(
	    "DSI", "Adding control points took: ",
	    std::chrono::duration_cast<std::chrono::seconds>(temp - duration)
		.count(),
	    " seconds");
	start = std::chrono::high_resolution_clock::now();
	int solver = 1;
	if (argc > 3) {
		solver = atoi(argv[3]);
	}
	dsi->solve_system(solver);
	// dsi->solve_system();
	duration = (std::chrono::high_resolution_clock::now() - start);
	Logger::out(
	    "DSI", "Solving system took: ",
	    std::chrono::duration_cast<std::chrono::seconds>(duration).count(),
	    " seconds");
	//// Tetrahedralize the GeoModel
	//// tetrahedralize(geomodel, NO_ID, false);
	// std::string file_name_out(ringmesh_test_data_path);
	GEO::vector<std::string> names;
	for (index_t m = 0; m < geomodel.nb_regions();
	     m++) {  // loop over regions
		const Region3D& cur_reg = geomodel.region(m);
		GEO::AttributesManager& reg_attr_mgr =
		    cur_reg.vertex_attribute_manager();

		reg_attr_mgr.list_attribute_names(names);

		for (const auto& name : names) {
			GEO::AttributeStore* attr_store =
			    reg_attr_mgr.find_attribute_store(name);
			GEO::ReadOnlyScalarAttributeAdapter cur_attr(
			    reg_attr_mgr, name);
			for (index_t vi = 0; vi < cur_reg.nb_vertices(); vi++) {
			}
		}
	}
	geomodel_save(geomodel, ringmesh_test_output_path + "geomodel_wdata" +
				    std::to_string(3) + "d." + "vtk");
}
