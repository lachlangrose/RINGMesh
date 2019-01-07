#ifndef MESHER_H
#define MESHER_H
#include <vector>
#include <ringmesh/geomodel/core/geomodel.h>
namespace RINGMesh {
class BoundingBoxMesher {
public:
	BoundingBoxMesher(GeoModel3D& geomodel);
	void mesh(std::vector<std::vector<double> > bounding_box,double max_volume);
private:
	GeoModel3D& geomodel_;
};
}



#endif /* MESHER_H */

