#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include <vector>

#include <ringmesh/geomodel/core/geomodel_mesh.h>
#include <Eigen/Dense>
#include <unordered_map>
namespace RINGMesh {

struct CellNeighbourIndex{
	CellNeighbourIndex(index_t cell1, index_t cell2) : _cell1(cell1),_cell2(cell2){}
	bool operator==(const CellNeighbourIndex& rhs) const;
	bool operator!=(const CellNeighbourIndex& rhs) const;
	index_t _cell1{ NO_ID };
	index_t _cell2{ NO_ID };
};
class Point3D {
public:
    Point3D(double x, double y, double z);
    inline double x() const;
    inline double y() const;
    inline double z() const;
private:
    double _x;
    double _y;
    double _z;
};
double Point3D::x() const {return _x;}
double Point3D::y() const {return _y;}
double Point3D::z() const {return _z;}
class InterfaceData : public Point3D {
public:
    InterfaceData(double x, double y, double z, double v);
    double _v;
    
private:
};
class PlanarData : public Point3D {
public:
    PlanarData(double x, double y, double z, Point3D v);
    Point3D _v;
};

class DSI {
public:
    DSI(const GeoModelMesh3D& mesh);
    void add_constant_gradient(double w = 0.1);
    void add_control_points(std::vector<InterfaceData> interfacedata);
    void add_gradient_control_points(std::vector<PlanarData> planardata);
    void build_interpolation_matrix();
    void solve_system();
private:
    void calculate_cell_gradient(index_t cell_global_index,
                                 Eigen::MatrixXd& grad,
                                 std::vector<index_t>& vert_idx);

    const GeoModelMesh3D& _mesh;
    std::vector<double> _A;
    std::vector<index_t> _row;
    std::vector<index_t> _col;
    std::vector<double> _B;

    Eigen::MatrixXd _I;
    index_t _c;
    index_t _nc;
};
}  // namespace RINGMesh
#endif /* INTERPOLATION_H */
