#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include <vector>

#include <unordered_map>
#include <ringmesh/geomodel/core/geomodel_mesh.h>
#include <Eigen/Dense>
namespace RINGMesh {
class Point3D {
public:
    Point3D(double x, double y, double z);

private:
    double _x;
    double _y;
    double _z;
};
class InterfaceData : Point3D {
public:
    InterfaceData(double x, double y, double z, double v);
    double _v;

private:
};
class PlanarData : Point3D {
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

private:
    void calculate_constant_gradient(index_t c1, index_t c2, std::vector<int>& idc, std::vector<double>& c);
    const GeoModelMesh3D& _mesh;
    std::vector<double> _A;
    std::vector<index_t> _row;
    std::vector<index_t> _col;
    std::unordered_map<index_t, index_t> _dinfo;
};
}
#endif /* INTERPOLATION_H */
