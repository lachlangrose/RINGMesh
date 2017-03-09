/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*!
 * @file Implementation of all GeoModelEntities classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geomodel/geomodel_mesh_entity.h>

#include <stack>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>
#include <ringmesh/geomodel/geomodel_validity.h>

namespace {
    using namespace RINGMesh ;

    /*!
     * @brief Checks that the geomodel vertex indices of @param E
     *       are in a valid range
     */
    bool check_range_model_vertex_ids( const GeoModelMeshEntity& E )
    {
        const GeoModelMeshVertices& geomodel_vertices = E.geomodel().mesh.vertices ;
        /// Check that the stored geomodel vertex indices are in a valid range
        for( index_t i = 0; i < E.nb_vertices(); ++i ) {
            if( geomodel_vertices.geomodel_vertex_id( E.gmme_id(), i ) == NO_ID
                && geomodel_vertices.geomodel_vertex_id( E.gmme_id(), i )
                    >= E.geomodel().mesh.vertices.nb() ) {
                Logger::warn( "GeoModelEntity" )
                    << "Invalid geomodel vertex index in " << E.gmme_id()
                    << std::endl ;
                return false ;
            }
        }
        return true ;
    }

    /*!
     * @brief Computes and returns the surface connected components
     * @details In Debug mode, the connected components are saved into 
     * an Attribute on surface facets.
     */
    index_t compute_nb_surface_connected_components( const Surface& M )
    {
        const index_t NO_COMPONENT = index_t( -1 ) ;
        GEO::Attribute< index_t > component( M.facet_attribute_manager(),
            "component" ) ;
        component.fill( NO_COMPONENT ) ;
        index_t nb_components = 0 ;
        for( index_t facet = 0; facet < M.nb_mesh_elements(); facet++ ) {
            if( component[facet] == NO_COMPONENT ) {
                std::stack< index_t > S ;
                S.push( facet ) ;
                component[facet] = nb_components ;
                do {
                    index_t cur_facet = S.top() ;
                    S.pop() ;
                    for( index_t edge = 0;
                        edge < M.nb_mesh_element_vertices( cur_facet ); edge++ ) {
                        index_t adj_facet = M.facet_adjacent_index( cur_facet,
                            edge ) ;
                        if( adj_facet != NO_ID
                            && component[adj_facet] == NO_COMPONENT ) {
                            S.push( adj_facet ) ;
                            component[adj_facet] = nb_components ;
                        }
                    }
                } while( !S.empty() ) ;
                nb_components++ ;
            }
        }
#ifndef RINGMESH_DEBUG
        component.destroy() ;
#endif
        return nb_components ;
    }

    /*!
     * @brief Computes and returns the region connected components
     * @details In Debug mode, the connected components are saved into 
     * an Attribute on region cells.
     */
    index_t compute_nb_volume_connected_components( const Region& M )
    {
        const index_t NO_COMPONENT = index_t( -1 ) ;
        GEO::Attribute< index_t > component( M.cell_attribute_manager(),
            "component" ) ;
        component.fill( NO_COMPONENT ) ;
        index_t nb_components = 0 ;
        for( index_t cell = 0; cell < M.nb_mesh_elements(); cell++ ) {
            if( component[cell] == NO_COMPONENT ) {
                std::stack< index_t > S ;
                S.push( cell ) ;
                component[cell] = nb_components ;
                do {
                    index_t cur_cell = S.top() ;
                    S.pop() ;
                    for( index_t facet = 0; facet < M.nb_cell_facets( cur_cell );
                        facet++ ) {
                        index_t adj_cell = M.cell_adjacent_index( cur_cell, facet ) ;
                        if( adj_cell != NO_ID
                            && component[adj_cell] == NO_COMPONENT ) {
                            S.push( adj_cell ) ;
                            component[adj_cell] = nb_components ;
                        }
                    }
                } while( !S.empty() ) ;
                nb_components++ ;
            }
        }
#ifndef RINGMESH_DEBUG
        component.destroy() ;
#endif
        return nb_components ;
    }

    /*!
     * @brief Count the number of times each vertex is in an edge or facet
     *
     * @param[in] gmme The GeoModelMeshEntity
     * @param[out] nb Resized to the number of vertices of the mesh.
     *      Number of times one vertex appear in an mesh_element collection of 
     *      the GeoModelMeshEntity edge or facet of the mesh.
     */
    void count_vertex_occurences(
        const GeoModelMeshEntity& E,
        std::vector< index_t >& nb )
    {
        nb.resize( E.nb_vertices(), 0 ) ;
        for( index_t mesh_element_index = 0;
            mesh_element_index < E.nb_mesh_elements(); ++mesh_element_index ) {
            for( index_t vertex = 0;
                vertex < E.nb_mesh_element_vertices( mesh_element_index );
                ++vertex ) {
                ++nb[E.mesh_element_vertex_index( mesh_element_index, vertex )] ;
            }
        }
    }

    index_t count_nb_isolated_vertices( const GeoModelMeshEntity& mesh )
    {
        std::vector< index_t > nb ;
        count_vertex_occurences( mesh, nb ) ;
        return static_cast< index_t >( std::count( nb.begin(), nb.end(), 0 ) ) ;
    }

    bool check_mesh_entity_vertices_are_different(
        std::vector< index_t >& vertices,
        std::vector< index_t >& vertices_global )
    {
        ringmesh_assert(
            std::count( vertices.begin(), vertices.end(), NO_ID ) == 0 ) ;
        ringmesh_assert(
            std::count( vertices_global.begin(), vertices_global.end(), NO_ID ) == 0 ) ;
        // 0 is the default value of the geomodel_vertex_id
        // If we have only 0 either this is a degenerate facets, but most certainly
        // geomodel vertex ids are not good
        ringmesh_assert(
            static_cast< index_t >( std::count( vertices_global.begin(), vertices_global.end(), 0 ) )
            != vertices_global.size() ) ;

        std::sort( vertices.begin(), vertices.end() ) ;
        std::sort( vertices_global.begin(), vertices_global.end() ) ;
        return std::unique( vertices.begin(), vertices.end() ) != vertices.end()
            || std::unique( vertices_global.begin(), vertices_global.end() )
                != vertices_global.end() ;
    }

    /*!
     * @brief Returns true if the surface facet is incident twice to the same vertex
     */
    bool facet_is_degenerate( const Surface& S, index_t f )
    {
        index_t nb_facet_vertices = S.nb_mesh_element_vertices( f ) ;
        std::vector< index_t > corners( nb_facet_vertices, NO_ID ) ;
        std::vector< index_t > corners_global( nb_facet_vertices, NO_ID ) ;
        index_t v = 0 ;
        const GeoModelMeshVertices& geomodel_vertices = S.geomodel().mesh.vertices ;
        for( index_t c = 0; c < S.nb_mesh_element_vertices( f ); ++c ) {
            index_t facet_vertex_index = S.mesh_element_vertex_index( f, c ) ;
            corners[v] = facet_vertex_index ;
            corners_global[v] = geomodel_vertices.geomodel_vertex_id( S.gmme_id(), f,
                v ) ;
            v++ ;
        }
        double area = S.mesh_element_size( f ) ;
        return check_mesh_entity_vertices_are_different( corners, corners_global )
            || area < S.geomodel().epsilon2() ;
    }

    /*!
     * @brief Returns true if the region cell is incident twice to the same vertex
     * or if the cell volume is negative or inferior to epsilon
     */
    bool cell_is_degenerate( const Region& region, index_t cell_index )
    {
        index_t nb_vertices_in_cell = region.nb_mesh_element_vertices( cell_index ) ;
        std::vector< index_t > vertices( nb_vertices_in_cell, NO_ID ) ;
        std::vector< index_t > vertices_global( nb_vertices_in_cell, NO_ID ) ;
        const GeoModelMeshVertices& geomodel_vertices =
            region.geomodel().mesh.vertices ;
        for( index_t v = 0; v < nb_vertices_in_cell; v++ ) {
            vertices[v] = region.mesh_element_vertex_index( cell_index, v ) ;
            vertices_global[v] = geomodel_vertices.geomodel_vertex_id(
                region.gmme_id(), cell_index, v ) ;
        }
        double volume = region.mesh_element_size( cell_index ) ;
        return check_mesh_entity_vertices_are_different( vertices, vertices_global )
            || volume < region.geomodel().epsilon3() ;
    }
}
/******************************************************************************/
namespace RINGMesh {
    bool GeoModelMeshEntity::is_inside_border( const GeoModelMeshEntity& rhs ) const
    {
        // Find out if this surface is twice in the in_boundary vector
        return std::count( in_boundary_.begin(), in_boundary_.end(), rhs.gmme_id() )
            > 1 ;
    }

    bool GeoModelMeshEntity::has_inside_border() const
    {
        for( index_t i = 0; i < nb_boundaries(); ++i ) {
            if( boundary( i ).is_inside_border( *this ) ) {
                return true ;
            }
        }
        return false ;
    }

    GeoModelMeshEntity::~GeoModelMeshEntity()
    {
        // Unbind attribute about vertex mapping
        GeoModel& modifiable_model = const_cast< GeoModel& >( geomodel() ) ;
        modifiable_model.mesh.vertices.unbind_geomodel_vertex_map( gmme_id() ) ;
#ifdef RINGMESH_DEBUG
        ringmesh_assert( mesh_ != NULL ) ;
        mesh_->print_mesh_bounded_attributes() ;
#endif
        delete mesh_ ;
    }

    bool GeoModelMeshEntity::are_geomodel_vertex_indices_valid() const
    {
        bool valid = true ;
        // For all vertices
        // Check that the global vertex has an index backward to 
        // the vertex of this entity
        const GeoModelMeshVertices& geomodel_vertices = geomodel().mesh.vertices ;
        for( index_t v = 0; v < nb_vertices(); ++v ) {
            index_t geomodel_v = geomodel_vertices.geomodel_vertex_id( gmme_id(),
                v ) ;

            if( geomodel_v == NO_ID ) {
                Logger::warn( "GeoModelEntity" ) << gmme_id() << " vertex " << v
                    << " is not mapped to the related global geomodel vertex indices."
                    << std::endl ;
                valid = false ;
            }

            std::vector< index_t > backward_vertices ;
            geomodel_vertices.mesh_entity_vertex_id( gmme_id(), geomodel_v,
                backward_vertices ) ;
            bool found_in_backward = false ;
            for( index_t bv = 0; bv < backward_vertices.size(); bv++ ) {
                if( backward_vertices[bv] == v ) {
                    found_in_backward = true ;
                }
            }
            if( !found_in_backward ) {
                Logger::warn( "GeoModelEntity" ) << "Error in mapping of "
                    << gmme_id() << " vertex " << v
                    << " to the related global geomodel vertex indices."
                    << std::endl ;
                valid = false ;
            }
        }
        return valid ;
    }

    bool GeoModelMeshEntity::is_index_valid() const
    {
        return index() < geomodel().nb_mesh_entities( type_name() ) ;
    }

    /*!
     * All entities in the boundary must have this in their
     *  in_boundary vector
     */
    bool GeoModelMeshEntity::is_boundary_connectivity_valid() const
    {
        const MeshEntityTypeManager& family = geomodel().entity_type_manager().mesh_entity_manager ;
        const MeshEntityType mesh_entity_type = type_name() ;
        const MeshEntityType& boundary_type = family.boundary_type( mesh_entity_type ) ;

        bool valid = true ;
        if( family.is_valid_type( boundary_type ) ) {
            for( index_t i = 0; i < nb_boundaries(); ++i ) {
                const GeoModelMeshEntity& E = boundary( i ) ;
                bool found = false ;
                index_t j = 0 ;
                while( !found && j < E.nb_in_boundary() ) {
                    if( E.in_boundary_gmme( j ) == gmme_id() ) {
                        found = true ;
                    }
                    j++ ;
                }
                if( !found ) {
                    Logger::warn( "GeoModelEntity" )
                        << "Inconsistency boundary-in_boundary between " << gmme_id()
                        << " and " << E.gmme_id() << std::endl ;
                    valid = false ;
                }
            }
        }
        return valid ;
    }
    /*! All entities must be at least in the boundary of another entity
     * and all entities in the in_boundary must have this entity in their
     * boundary vector
     */
    bool GeoModelMeshEntity::is_in_boundary_connectivity_valid() const
    {
        const MeshEntityTypeManager& family = geomodel().entity_type_manager().mesh_entity_manager ;
        const MeshEntityType mesh_entity_type = type_name() ;
        const MeshEntityType& in_boundary_type = family.in_boundary_type( mesh_entity_type ) ;

        bool valid = true ;
        if( family.is_valid_type( in_boundary_type ) ) {
            if( nb_in_boundary() == 0 ) {
                Logger::warn( "GeoModelEntity" ) << gmme_id()
                    << " is in the boundary of no entity " << std::endl ;
                valid = false ;
            }
            for( index_t i = 0; i < nb_in_boundary(); ++i ) {
                const GeoModelMeshEntity& E = in_boundary( i ) ;
                bool found = false ;
                index_t j = 0 ;
                while( !found && j < E.nb_boundaries() ) {
                    if( E.boundary_gmme( j ) == gmme_id() ) {
                        found = true ;
                    }
                    j++ ;
                }
                if( !found ) {
                    Logger::warn( "GeoModelEntity" )
                        << "Inconsistency in_boundary-boundary between " << gmme_id()
                        << " and " << E.gmme_id() << std::endl ;
                    valid = false ;
                }
            }
        }
        return valid ;
    }
    /*!
     *  If the parent type is defined for this EntityType, 
     *  and if the geomodel has entities of that type, the entity must have a parent.
     * @todo Remove the second if condition ?
     */
    bool GeoModelMeshEntity::is_parent_connectivity_valid() const
    {
        bool valid = true ;

        const RelationshipManager& family = geomodel().entity_type_manager().relationship_manager ;
        const MeshEntityType mesh_entity_type = type_name() ;

        const std::vector< GeologicalEntityType > parent_types = family.parent_types(
            mesh_entity_type ) ;
        for( index_t p_itr = 0; p_itr < parent_types.size(); ++p_itr ) {
            const GeologicalEntityType& parent_type = parent_types[p_itr] ;
            index_t nb_parent_entities_in_geomodel =
                geomodel_.nb_geological_entities( parent_type ) ;
            if( nb_parent_entities_in_geomodel == 0 ) {
                continue ;
            } else {
                // There must be one and only one parent of that type in this entity
                // And this parent msut have this entity in its children
                index_t nb_found_parents = 0 ;
                for( index_t i = 0; i < nb_parents(); ++i ) {
                    const GeoModelGeologicalEntity& GMGE = parent( i ) ;
                    if( GMGE.type_name() == parent_type ) {
                        nb_found_parents++ ;

                        // The parent must have this entity in its children
                        bool found = false ;
                        index_t j = 0 ;
                        while( !found && j < GMGE.nb_children() ) {
                            if( GMGE.child_gmme( j ) == gmme_id() ) {
                                found = true ;
                            }
                            j++ ;
                        }
                        if( !found ) {
                            Logger::warn( "GeoModelEntity" )
                                << "Inconsistency parent-child between " << gmme_id()
                                << " and " << GMGE.gmge_id() << std::endl ;
                            valid = false ;
                        }
                    }
                }
                if( nb_found_parents != 1 ) {
                    Logger::warn( "GeoModelEntity" ) << gmme_id() << " has "
                        << nb_found_parents << " geological parent entity of type "
                        << parent_type << std::endl ;
                    valid = false ;
                }
            }
        }
        return valid ;
    }
    /*!  Check that required information for the TYPE is defined
     *    and that reverse information is stored by the corresponding
     *    entities
     */
    bool GeoModelMeshEntity::is_connectivity_valid() const
    {
        return is_boundary_connectivity_valid()
            && is_in_boundary_connectivity_valid() && is_parent_connectivity_valid() ;
    }


    bool GeoModelMeshEntity::is_identification_valid() const
    {
        bool is_valid = true ;
        if( !gmme_id().is_defined() ) {
            Logger::err( "GeoModeMeshlEntity" ) << " Entity associated to geomodel "
                << geomodel().name() << "has no type and/or no index " << std::endl ;
            is_valid = false ;
            // No further checks are possible - This really should not happen
            ringmesh_assert_not_reached;
        }
        if( !is_index_valid() ) {
            Logger::warn( "GeoModelMeshEntity" ) << " Entity index " << gmme_id()
                << " is not valid. " << std::endl ;
            // This really should not happen
            is_valid = false ;
            ringmesh_assert_not_reached;
        }
        return is_valid ;
    }

    /*!
     * @return Assert that the parent exists and returns it.
     */
    const GeoModelGeologicalEntity& GeoModelMeshEntity::parent(
        index_t parent_index ) const
    {
        gmge_t parent = parent_gmge( parent_index ) ;
        ringmesh_assert( parent.is_defined() ) ;
        return geomodel().geological_entity( parent ) ;
    }
    const GeoModelGeologicalEntity& GeoModelMeshEntity::parent(
        const GeologicalEntityType& parent_type ) const
    {
        gmge_t id = parent_gmge( parent_type ) ;
        ringmesh_assert( id.is_defined() ) ;
        return geomodel().geological_entity( id ) ;
    }
    const gmge_t& GeoModelMeshEntity::parent_gmge(
        const GeologicalEntityType& parent_type_name ) const
    {
        for( index_t i = 0; i < nb_parents(); ++i ) {
            if( parents_[i].type() == parent_type_name ) {
                return parent_gmge( i ) ;
            }
        }
        ringmesh_assert_not_reached ;
    }

    const GeoModelMeshEntity& GeoModelMeshEntity::boundary( index_t x ) const
    {
        return geomodel().mesh_entity( boundary_gmme( x ) ) ;
    }

    const GeoModelMeshEntity& GeoModelMeshEntity::in_boundary( index_t x ) const
    {
        return geomodel().mesh_entity( in_boundary_gmme( x ) ) ;
    }

    /**************************************************************/
    /*!
     * @brief Checks if this entity define the geomodel external boundary
     * @details Test if the entity is in the Surfaces defining the universe
     */
    bool Corner::is_on_voi() const
    {
        // True if one of the incident surface define the universe
        for( index_t i = 0; i < nb_in_boundary(); ++i ) {
            if( in_boundary( i ).is_on_voi() ) {
                return true ;
            }
        }
        return false ;

    }
    /*!
     * @brief Check that the Corner mesh is a unique point
     */
    bool Corner::is_mesh_valid() const
    {
        bool valid = true ;
        if( nb_vertices() != 1 ) {
            Logger::err( "GeoModelEntity" ) << "Corner " << index() << " mesh has "
                << mesh0d_->nb_vertices() << " vertices " << std::endl ;
            valid = false ;
        }

        // The default point is (0., 0., 0.) and there might be a valid
        // Corner at this position.
        /*if( mesh_.vertices.point( 0 ) == vec3() ) {
         Logger::warn( "GeoModelEntity" )
         << "Corner " << index()
         << " point is default " << std::endl ;
         valid = false ;
         }*/
        return valid ;
    }

    /***************************************************************/

    /*!
     * @brief Construct a Line
     *
     * @param[in] geomodel The parent geomodel
     * @param[in] id The index of the line in the lines_ vector of the parent geomodel
     */
    Line::Line( const GeoModel& geomodel, index_t id, const MeshType type )
        : GeoModelMeshEntity( geomodel, id )
    {
        update_mesh_storage_type( Mesh1D::create_mesh( type ) ) ;
    }

    /*!
     * @brief Check that the mesh of the Line is valid
     * @details Check that 
     *  - the GEO::Mesh has more than 1 vertex - more than 1 edge - no facets - no cells.
     *  - global indices of vertices in the geomodel are in a valid range 
     *  - each vertex is in 2 edges except extremities that are in 1 edge
     * 
     * Does not check:
     *  - Self-intersection - I suppose there are no segment - segment intersection (JP)
     *  - Duplicated edge - most probably ruled out with the duplicated vertex test (JP)
     *  - Duplicated vertex (verified at GeoModel level)
     */
    bool Line::is_mesh_valid() const
    {
        bool valid = true ;

        // Check that the GEO::Mesh has the expected entities
        if( nb_vertices() < 2 ) {
            Logger::err( "GeoModelEntity" ) << "Line " << index() << " has "
                << mesh1d_->nb_vertices() << " vertices " << std::endl ;
            valid = false ;
        }
        if( mesh1d_->nb_edges() == 0 ) {
            Logger::err( "GeoModelEntity" ) << "Line " << index() << " mesh has "
                << mesh1d_->nb_edges() << " edges " << std::endl ;
            valid = false ;
        }

        // Model indices must be valid
        valid = check_range_model_vertex_ids( *this ) && valid ;

        if( nb_vertices() > 1 ) {
            // Count the number of edges in which each vertex is
            std::vector< index_t > nb ;
            count_vertex_occurences( *this, nb ) ;
            index_t nb0 = 0 ;
            index_t nb1 = 0 ;
            index_t nb2 = 0 ;
            for( index_t i = 0; i < nb.size(); ++i ) {
                if( nb[i] == 0 )
                    ++nb0 ;
                else if( nb[i] == 1 )
                    ++nb1 ;
                else if( nb[i] == 2 ) ++nb2 ;
            }

            // Vertices at extremitites must be in only one edge
            if( nb.front() != 1 || nb.back() != 1 ) {
                Logger::err( "GeoModelEntity" ) << "Invalid extremity points in "
                    << gmme_id() << std::endl ;
                valid = false ;
            }
            // No isolated vertices are allowed
            if( nb0 > 0 ) {
                Logger::warn( "GeoModelEntity" ) << nb0 << " isolated vertices in "
                    << gmme_id() << std::endl ;
                valid = false ;
            }
            // Only the two extremities are in only 1 edge 
            // One connected component condition
            if( nb1 != 2 ) {
                Logger::warn( "GeoModelEntity" )
                    << "More than one connected component for " << gmme_id()
                    << std::endl ;
                valid = false ;
            }
            // All the others must be in 2 edges and 2 edges only
            // Manifold condition
            if( nb2 != nb.size() - 2 ) {
                Logger::warn( "GeoModelEntity" ) << "Non-manifold entity" << gmme_id()
                    << std::endl ;
                valid = false ;
            }
        }

        // No zero edge length
        index_t nb_degenerated = 0 ;
        for( index_t e = 0; e < nb_mesh_elements(); ++e ) {
            double l = length(
                mesh_element_vertex( e, 1 ) - mesh_element_vertex( e, 0 ) ) ;
            if( l < geomodel().epsilon() ) {
                nb_degenerated++ ;
            }
        }
        if( nb_degenerated > 0 ) {
            Logger::warn( "GeoModelEntity" ) << nb_degenerated
                << " degenerated edges in " << gmme_id() << std::endl ;
            valid = false ;
        }
        if( !is_first_corner_first_vertex() ) {
            Logger::warn( "GeoModelEntity" ) << "First and last vertex of Line"
                << index() << " do not match first and second Corner respectively"
                << std::endl ;
            valid = false ;
        }
        return valid ;
    }

    bool Line::is_connectivity_valid() const
    {
        bool line_valid = GeoModelMeshEntity::is_connectivity_valid() ;

        // A Line must have 2 corners - they are identical if the Line is closed
        if( nb_boundaries() != 2 ) {
            Logger::warn( "GeoModelEntity" ) << gmme_id()
                << " does not have 2 corners" << std::endl ;
            line_valid = false ;
        }
        return line_valid ;

    }

    bool Line::is_first_corner_first_vertex() const
    {
        if( nb_boundaries() != 2 || nb_vertices() < 2 ) {
            return false ;
        } else {
            // Geometric comparison - not great at all
            return boundary( 0 ).vertex( 0 ) == vertex( 0 )
                && boundary( 1 ).vertex( 0 ) == vertex( nb_vertices() - 1 ) ;
        }
    }

    bool Line::is_on_voi() const
    {
        // True if one of the incident surface define the universe
        for( index_t i = 0; i < nb_in_boundary(); ++i ) {
            if( in_boundary( i ).is_on_voi() ) {
                return true ;
            }
        }
        return false ;
    }

    /********************************************************************/

    /*!
     * @brief Check that the mesh of the Surface is valid
     * @details Check that
     *  - the GEO::Mesh has more than 2 vertices, at least 1 facet, no cells.
     *  - global indices of vertices in the geomodel are in a valid range
     *  - no degenerate facet 
     *  - one connected component 
     *
     *  Some tests are not performed here but globally on the GeoModel
     *  - intersection of facets 
     *  - non-manifold edges 
     *  - duplicated vertices are on a boundary Line ending in the Surface 
     * 
     *
     *  Some tests are not performed     
     *  - non-manifold points
     *  - surface orientability
     *  - planarity of polygonal facets 
     *
     * @todo Check that there is no duplicated facet 
     */
    bool Surface::is_mesh_valid() const
    {
        bool valid = true ;
        // Check that the GEO::Mesh has the expected entities
        // at least 3 vertices and one facet.
        if( nb_vertices() < 3 ) {
            Logger::warn( "GeoModelEntity" ) << gmme_id()
                << " has less than 3 vertices " << std::endl ;
            valid = false ;
        }
        if( mesh2d_->nb_facets() == 0 ) {
            Logger::warn( "GeoModelEntity" ) << gmme_id() << " has no facets "
                << std::endl ;
            valid = false ;
        }

        // No isolated vertices
        index_t nb_isolated_vertices = count_nb_isolated_vertices( *this ) ;
        if( nb_isolated_vertices > 0 ) {
            Logger::warn( "GeoModelEntity" ) << gmme_id() << " mesh has "
                << nb_isolated_vertices << " isolated vertices " << std::endl ;
            valid = false ;
        }

        // No zero area facet
        // No facet incident to the same vertex check local and global indices
        index_t nb_degenerate = 0 ;
        for( index_t f = 0; f < mesh2d_->nb_facets(); f++ ) {
            if( facet_is_degenerate( *this, f ) ) {
                nb_degenerate++ ;
            }
        }
        if( nb_degenerate != 0 ) {
            Logger::warn( "GeoModelEntity" ) << gmme_id() << " mesh has "
                << nb_degenerate << " degenerate facets " << std::endl ;
            valid = false ;
        }

        // No duplicated facet
        GEO::vector< index_t > colocated ;
        // GEO::mesh_detect_duplicated_facets( mesh_, colocated ) ; // not implemented yet 
        index_t nb_duplicated_f = 0 ;
        for( index_t f = 0; f < colocated.size(); ++f ) {
            if( colocated[f] != f ) {
                nb_duplicated_f++ ;
            }
        }
        if( nb_duplicated_f > 0 ) {
            Logger::warn( "GeoModelEntity" ) << gmme_id() << " mesh has "
                << nb_duplicated_f << " duplicated facets " << std::endl ;
            valid = false ;
        }

        // One connected component  
        index_t cc = compute_nb_surface_connected_components( *this ) ;
        if( cc != 1 ) {
            Logger::warn( "GeoModelEntity" ) << gmme_id() << " mesh has " << cc
                << " connected components " << std::endl ;
            valid = false ;
        }
        return valid ;
    }

    bool Surface::is_on_voi() const
    {
        for( index_t i = 0; i < geomodel().universe().nb_boundaries(); ++i ) {
            if( geomodel().universe().boundary_gmme( i ) == gmme_id() ) {
                return true ;
            }
        }
        return false ;
    }

    void Surface::next_on_border(
        index_t f,
        index_t e,
        index_t& next_f,
        index_t& next_e ) const
    {
        ringmesh_assert( e < nb_mesh_element_vertices( f ) ) ;
        ringmesh_assert( is_on_border( f, e ) ) ;

        // Global indices in the surfaces
        index_t next_v_id = mesh_element_vertex_index( f,
            next_facet_vertex_index( f, e ) ) ;

        // Get the facets around the shared vertex (next_v_id) that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)
        std::vector< index_t > facets_around_next_v_id ;
        index_t nb_around = facets_around_vertex( next_v_id, facets_around_next_v_id,
            true, f ) ;
        ringmesh_assert( nb_around == 1 || nb_around == 2 ) ;

        next_f = facets_around_next_v_id[0] ;

        if( nb_around == 2 ) {
            if( next_f == f ) {
                next_f = facets_around_next_v_id[1] ;
            }
            ringmesh_assert( next_f != NO_ID ) ;
            ringmesh_assert( is_on_border( next_f ) ) ;

            // Local index of next vertex in the next facet
            next_e = vertex_index_in_facet( next_f, next_v_id ) ;
            ringmesh_assert( is_on_border( next_f, next_e ) ) ;
        } else if( nb_around == 1 ) {
            // next_v_id must be in two border edges of facet f
            next_e = vertex_index_in_facet( next_f, next_v_id ) ;
            ringmesh_assert( is_on_border( next_f, next_e ) ) ;
        }
    }

    void Surface::prev_on_border(
        index_t f,
        index_t e,
        index_t& prev_f,
        index_t& prev_e ) const
    {
        ringmesh_assert( e < nb_mesh_element_vertices( f ) ) ;
        ringmesh_assert( is_on_border( f, e ) ) ;

        // Global indices in the surfaces
        index_t v_id = mesh_element_vertex_index( f, e ) ;

        // Get the facets around the shared vertex (v_id) that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)
        std::vector< index_t > facets_around_v_id ;
        index_t nb_around = facets_around_vertex( v_id, facets_around_v_id, true,
            f ) ;
        ringmesh_assert( nb_around == 1 || nb_around == 2 ) ;

        prev_f = facets_around_v_id[0] ;

        if( nb_around == 2 ) {
            if( prev_f == f ) {
                prev_f = facets_around_v_id[1] ;
            }
            ringmesh_assert( prev_f != NO_ID ) ;
            ringmesh_assert( is_on_border( prev_f ) ) ;

            // Local index of given vertex in the prev facet
            index_t v_in_prev_f = vertex_index_in_facet( prev_f, v_id ) ;
            // Local index of previous vertex in the prev facet
            prev_e = prev_facet_vertex_index( prev_f, v_in_prev_f ) ;
            ringmesh_assert( is_on_border( prev_f, prev_e ) ) ;
        } else if( nb_around == 1 ) {
            // v_id must be in two border edges of facet f
            index_t v_in_next_facet = vertex_index_in_facet( prev_f, v_id ) ;
            prev_e = prev_facet_vertex_index( prev_f, v_in_next_facet ) ;
            ringmesh_assert( is_on_border( prev_f, prev_e ) ) ;
        }
    }

    /*!
     * @brief Get the first facet of the surface that has an edge linking the two vertices (ids in the surface)
     *
     * @param[in] in0 Index of the first vertex in the surface
     * @param[in] in1 Index of the second vertex in the surface
     * @return NO_ID or the index of the facet
     */
    index_t Surface::facet_from_surface_vertex_ids( index_t in0, index_t in1 ) const
    {
        ringmesh_assert(
            in0 < nb_vertices() && in1 <nb_vertices() ) ;

        // Another possible, probably faster, algorithm is to check if the 2 indices
        // are neighbors in facets_ and check that they are in the same facet

        // Check if the edge is in one of the facet
        for( index_t f = 0; f < nb_mesh_elements(); ++f ) {
            bool found = false ;
            index_t prev = mesh_element_vertex_index( f,
                nb_mesh_element_vertices( f ) - 1 ) ;
            for( index_t v = 0; v < nb_mesh_element_vertices( f ); ++v ) {
                index_t p = mesh_element_vertex_index( f, v ) ;
                if( ( prev == in0 && p == in1 ) || ( prev == in1 && p == in0 ) ) {
                    found = true ;
                    break ;
                }
                prev = p ;
            }
            if( found ) {
                return f ;
            }
        }
        return NO_ID ;
    }

    /*!
     * @brief Comparator of two vec3
     */
    struct comp_vec3bis {
        bool operator()( const vec3& l, const vec3& r ) const
        {
            if( l.x != r.x ) {
                return l.x < r.x ;
            }
            if( l.y != r.y ) {
                return l.y < r.y ;
            }
            return l.z < r.z ;
        }
    } ;

    index_t Surface::facets_around_vertex(
        index_t surf_vertex_id,
        std::vector< index_t >& result,
        bool border_only,
        index_t f0 ) const
    {
        result.clear() ;

        index_t f = 0 ;
        while( f0 == NO_ID && f < nb_mesh_elements() ) {
            for( index_t lv = 0; lv < nb_mesh_element_vertices( f ); lv++ ) {
                if( mesh_element_vertex_index( f, lv ) == surf_vertex_id ) {
                    f0 = f ;
                    break ;
                }
            }
            f++ ;
        }

        ringmesh_assert( f0 != NO_ID ) ;

        // Flag the visited facets
        std::vector< index_t > visited ;
        visited.reserve( 10 ) ;

        // Stack of the adjacent facets
        std::stack< index_t > S ;
        S.push( f0 ) ;
        visited.push_back( f0 ) ;

        do {
            index_t f = S.top() ;
            S.pop() ;

            for( index_t v = 0; v < nb_mesh_element_vertices( f ); ++v ) {
                if( mesh_element_vertex_index( f, v ) == surf_vertex_id ) {
                    index_t adj_P = facet_adjacent_index( f, v ) ;
                    index_t prev = prev_facet_vertex_index( f, v ) ;
                    index_t adj_prev = facet_adjacent_index( f, prev ) ;

                    if( adj_P != NO_ID ) {
                        // The edge starting at P is not on the boundary
                        if( !contains( visited, adj_P ) ) {
                            S.push( adj_P ) ;
                            visited.push_back( adj_P ) ;
                        }
                    }
                    if( adj_prev != NO_ID ) {
                        // The edge ending at P is not on the boundary
                        if( !contains( visited, adj_prev ) ) {
                            S.push( adj_prev ) ;
                            visited.push_back( adj_prev ) ;
                        }
                    }

                    if( border_only ) {
                        if( adj_P == NO_ID || adj_prev == NO_ID ) {
                            result.push_back( f ) ;
                        }
                    } else {
                        result.push_back( f ) ;
                    }

                    // We are done with this facet
                    break ;
                }
            }
        } while( !S.empty() ) ;

        return static_cast< index_t >( result.size() ) ;
    }

    /*!
     * @brief Compute closest vertex in a facet to a point
     * @param[in] f Facet index
     * @param[in] v Coordinates of the point to which distance is measured
     * @return Index of the vertex of @param f closest to @param v
     */
    index_t Surface::closest_vertex_in_facet( index_t f, const vec3& v ) const
    {
        index_t result = 0 ;
        double dist = DBL_MAX ;
        for( index_t p = 0; p < nb_mesh_element_vertices( f ); p++ ) {
            double distance = length2( v - mesh_element_vertex( f, p ) ) ;
            if( dist > distance ) {
                dist = distance ;
                result = p ;
            }
        }
        return result ;
    }

    /********************************************************************/

    bool Region::is_on_voi() const
    {
        return false ;
    }

    bool Region::is_connectivity_valid() const
    {
        if( nb_boundaries() != sides_.size() ) {
            Logger::err( "GeoModelEntity" ) << gmme_id()
                << " boundary sides are invalid " << std::endl ;
            ringmesh_assert_not_reached ;
            return false ;
        }
        bool region_valid = GeoModelMeshEntity::is_connectivity_valid() ;
        if( nb_boundaries() == 0 ) {
            Logger::warn( "GeoModelEntity" ) << gmme_id() << " has no boundaries "
                << std::endl ;
            region_valid = false ;
        }
        return region_valid ;
    }

    bool Region::is_mesh_valid() const
    {
        if( !is_meshed() ) {
            return true ;
        } else {
            bool valid = true ;
            // Check that the GEO::Mesh has the expected entities
            // at least 4 vertices and one cell.
            if( mesh3d_->nb_vertices() < 4 ) {
                Logger::warn( "GeoModelEntity" ) << gmme_id()
                    << " has less than 4 vertices " << std::endl ;
                valid = false ;
            }

            // No isolated vertices
            index_t nb_isolated_vertices = count_nb_isolated_vertices( *this ) ;
            if( nb_isolated_vertices > 0 ) {
                Logger::warn( "GeoModelEntity" ) << gmme_id() << " mesh has "
                    << nb_isolated_vertices << " isolated vertices " << std::endl ;
                valid = false ;
            }

            // No cell with negative volume
            // No cell incident to the same vertex check local and global indices
            index_t nb_degenerate = 0 ;
            for( index_t c = 0; c < mesh3d_->nb_cells(); c++ ) {
                if( cell_is_degenerate( *this, c ) ) {
                    nb_degenerate++ ;
                }
            }
            if( nb_degenerate != 0 ) {
                Logger::warn( "GeoModelEntity" ) << gmme_id() << " mesh has "
                    << nb_degenerate << " degenerate cells " << std::endl ;
                valid = false ;
            }

            // One connected component
            index_t cc = compute_nb_volume_connected_components( *this ) ;
            if( cc != 1 ) {
                Logger::warn( "GeoModelEntity" ) << gmme_id() << " mesh has " << cc
                    << " connected components " << std::endl ;
                valid = false ;
            }
            return valid ;
        }
    }

    void Region::compute_region_volumes_per_cell_type(
        double& tet_volume,
        double& pyramid_volume,
        double& prism_volume,
        double& hex_volume,
        double& poly_volume ) const
    {
        for( index_t c = 0; c < nb_mesh_elements(); c++ ) {
            index_t nb_vertices = nb_mesh_element_vertices( c ) ;
            double volume = mesh3d_->cell_volume( c ) ;
            switch( nb_vertices ) {
                case 4:
                    tet_volume += volume ;
                    break ;
                case 5:
                    pyramid_volume += volume ;
                    break ;
                case 6:
                    prism_volume += volume ;
                    break ;
                case 8:
                    hex_volume += volume ;
                    break ;
                default:
                    poly_volume += volume ;
                    break ;
            }
        }
    }

    index_t Region::cells_around_vertex(
        index_t vertex_id,
        std::vector< index_t >& result,
        index_t cell_hint ) const
    {
        result.resize( 0 ) ;

        if( cell_hint == NO_ID ) {
            return 0 ;
        }

        // Flag the visited cells
        std::vector< index_t > visited ;
        visited.reserve( 10 ) ;

        // Stack of the adjacent cells
        std::stack< index_t > S ;
        S.push( cell_hint ) ;
        visited.push_back( cell_hint ) ;

        do {
            index_t c = S.top() ;
            S.pop() ;

            bool cell_includes_vertex = false ;
            for( index_t v = 0; v < nb_mesh_element_vertices( c ); v++ ) {
                if( mesh_element_vertex_index( c, v ) == vertex_id ) {
                    result.push_back( c ) ;
                    cell_includes_vertex = true ;
                    break ;
                }
            }
            if( !cell_includes_vertex ) {
                continue ;
            }

            for( index_t f = 0; f < nb_cell_facets( c ); f++ ) {
                for( index_t v = 0; v < nb_cell_facet_vertices( c, f ); v++ ) {
                    index_t vertex = cell_facet_vertex_index( c, f, v ) ;
                    if( vertex == vertex_id ) {
                        index_t adj_P = cell_adjacent_index( c, f ) ;

                        if( adj_P != NO_ID ) {
                            if( !contains( visited, adj_P ) ) {
                                S.push( adj_P ) ;
                                visited.push_back( adj_P ) ;
                            }
                        }
                        break ;
                    }
                }
            }
        } while( !S.empty() ) ;

        return static_cast< index_t >( result.size() ) ;
    }

    void GeoModelMeshEntityAccess::change_mesh_data_structure( const MeshType type )
    {
        if( MeshEntityTypeManager::is_corner( gmme_.type_name() ) ) {
            Corner& corner = dynamic_cast< Corner& >( gmme_ ) ;
            Mesh0D* old_mesh = corner.mesh0d_ ;
            if( old_mesh->type_name() == type ) {
                return ;
            }
            corner.update_mesh_storage_type( Mesh0D::create_mesh( type ) ) ;
            Mesh0DBuilder_var builder = Mesh0DBuilder::create_builder(
                *corner.mesh0d_ ) ;
            builder->copy( *old_mesh, true ) ;
            delete old_mesh ;
        } else if( MeshEntityTypeManager::is_line( gmme_.type_name() ) ) {
            Line& line = dynamic_cast< Line& >( gmme_ ) ;
            Mesh1D* old_mesh = line.mesh1d_ ;
            if( old_mesh->type_name() == type ) {
                return ;
            }
            line.update_mesh_storage_type( Mesh1D::create_mesh( type ) ) ;
            Mesh1DBuilder_var builder = Mesh1DBuilder::create_builder(
                *line.mesh1d_ ) ;
            builder->copy( *old_mesh, true ) ;
            delete old_mesh ;
        } else if( MeshEntityTypeManager::is_surface( gmme_.type_name() ) ) {
            Surface& surface = dynamic_cast< Surface& >( gmme_ ) ;
            Mesh2D* old_mesh = surface.mesh2d_ ;
            if( old_mesh->type_name() == type ) {
                return ;
            }
            surface.update_mesh_storage_type( Mesh2D::create_mesh( type ) ) ;
            Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
                *surface.mesh2d_ ) ;
            builder->copy( *old_mesh, true ) ;
            delete old_mesh ;
        } else if( MeshEntityTypeManager::is_region( gmme_.type_name() ) ) {
            Region& region = dynamic_cast< Region& >( gmme_ ) ;
            Mesh3D* old_mesh = region.mesh3d_ ;
            if( old_mesh->type_name() == type ) {
                return ;
            }
            region.update_mesh_storage_type( Mesh3D::create_mesh( type ) ) ;
            Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
                *region.mesh3d_ ) ;
            builder->copy( *old_mesh, true ) ;
            delete old_mesh ;
        } else {
            ringmesh_assert_not_reached ;
        }
    }
}
