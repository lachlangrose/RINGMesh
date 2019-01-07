#include <ringmesh/mesh/mesh_builder.h>
#include <geogram/mesh/mesh.h>

#include <ringmesh/ringmesh_tests_config.h>
#include <mesher.h>

#include <geogram/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder_from_mesh.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/mesh_set.h>
using namespace RINGMesh;

BoundingBoxMesher::BoundingBoxMesher(GeoModel3D& geomodel) : geomodel_(geomodel){}

void BoundingBoxMesher::mesh(std::vector<std::vector<double> > bounding_box,double max_volume) {
	GeoModelBuilder3D builder( geomodel_ );
	vec3 v0( bounding_box[0][0], bounding_box[1][0], bounding_box[2][0] );
    	vec3 v1( bounding_box[0][1], bounding_box[1][0], bounding_box[2][0] );
    	vec3 v2( bounding_box[0][1], bounding_box[1][1], bounding_box[2][0] );
    	vec3 v3( bounding_box[0][0], bounding_box[1][1], bounding_box[2][0] );
    	vec3 v4( bounding_box[0][0], bounding_box[1][0], bounding_box[2][1] );
    	vec3 v5( bounding_box[0][1], bounding_box[1][0], bounding_box[2][1]);
    	vec3 v6( bounding_box[0][1], bounding_box[1][1], bounding_box[2][1]);
    	vec3 v7( bounding_box[0][0], bounding_box[1][1], bounding_box[2][1]);
    	std::vector< index_t > triangles;
    	triangles.push_back( 0 );
    	triangles.push_back( 1 );
    	triangles.push_back( 2 );
    	triangles.push_back( 0 );
    	triangles.push_back( 2 );
    	triangles.push_back( 3 );

    	std::vector< index_t > surface_facet_ptr;
    	surface_facet_ptr.push_back( 0 );
    	surface_facet_ptr.push_back( 3 );
    	surface_facet_ptr.push_back( 6 );

    	std::vector< vec3 > vertices( 4 );
    	gmme_id id;

    	id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    	vertices[0] = v0;
    	vertices[1] = v1;
    	vertices[2] = v2;
    	vertices[3] = v3;
    	builder.geometry.set_surface_geometry(
    	    id.index(), vertices, triangles, surface_facet_ptr );

    	id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    	vertices[0] = v1;
    	vertices[1] = v5;
    	vertices[2] = v6;
    	vertices[3] = v2;
    	builder.geometry.set_surface_geometry(
    	    id.index(), vertices, triangles, surface_facet_ptr );

    	id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    	vertices[0] = v4;
    	vertices[1] = v5;
    	vertices[2] = v6;
    	vertices[3] = v7;
    	builder.geometry.set_surface_geometry(
    	    id.index(), vertices, triangles, surface_facet_ptr );

    	id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    	vertices[0] = v0;
    	vertices[1] = v3;
    	vertices[2] = v7;
    	vertices[3] = v4;
    	builder.geometry.set_surface_geometry(
    	    id.index(), vertices, triangles, surface_facet_ptr );

    	id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    	vertices[0] = v3;
    	vertices[1] = v2;
    	vertices[2] = v6;
    	vertices[3] = v7;
    	builder.geometry.set_surface_geometry(
    	    id.index(), vertices, triangles, surface_facet_ptr );

    	id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    	vertices[0] = v0;
    	vertices[1] = v1;
    	vertices[2] = v5;
    	vertices[3] = v4;
    	builder.geometry.set_surface_geometry(
        id.index(), vertices, triangles, surface_facet_ptr );
        builder.build_lines_and_corners_from_surfaces();
	builder.build_regions_from_lines_and_surfaces();
        builder.end_geomodel();
     	
	//// Tetrahedralize the GeoModel
        tetrahedralize( geomodel_, NO_ID, true ,max_volume);

}
int main() {
   GeoModel3D geomodel;
   std::vector<std::vector<double> > bb;
   bb.push_back({0., 2.});
   bb.push_back({0., 2.});
   bb.push_back({0., 2.});
   BoundingBoxMesher mesher(geomodel);

   mesher.mesh(bb,0.2);
       geomodel_save(geomodel, ringmesh_test_output_path + "geomodel_wdata" +
                                std::to_string(3) + "d." + "vtk");
}

