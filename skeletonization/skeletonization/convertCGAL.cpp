//
//  convertCGAL.cpp
//  skeletonization
//
//  Pre-condition: Object's translate must be set to 0,0,0
//
//  Created by Gary Tse on 5/4/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef convertCGAL_cpp
#define convertCGAL_cpp

//#include "convertCGAL.hpp"

#include <cstdio>

#include <CGAL/Point_3.h>
#include <CGAL/Cartesian.h>
//#include <CGAL/HalfedgeDS_face_max_base_with_id.h>
#include <CGAL/HalfedgeDS_face_base.h>
//#include <CGAL/HalfedgeDS_vertex_max_base_with_id.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
//#include <CGAL/Mean_curvature_flow_skeletonization.h>

#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>

template <class Refs>
struct My_face : public CGAL::HalfedgeDS_face_base<Refs> {
    unsigned id;
};

template <class Refs, class Point>
struct My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> {
    My_vertex() {}
    My_vertex(const Point& pt) : CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>(pt) {}
    unsigned id;
};

// An items type using my face and vertex.
struct My_items : public CGAL::Polyhedron_items_3 {
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef My_face<Refs> Face;
    };
    
    template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point;
        typedef My_vertex<Refs, Point> Vertex;
    };
    
};

typedef CGAL::Cartesian<double> Kernel;
//typedef CGAL::Polyhedron_3<Kernel, My_items> Polyhedron;                          //original using custom Face & Vertex IDs
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

// A modifier converting a Maya mesh to CGAL Polyhedron_3
template <class HDS>
class Mesh_to_polyhedron : public CGAL::Modifier_base<HDS> {
public:
    Mesh_to_polyhedron(MFnMesh &mesh) : m_mesh(mesh.object()) {}
    void operator()(HDS& hds) {
        // get mesh data
        MFloatPointArray pts;
        m_mesh.getPoints(pts);
        MIntArray tcounts, tvers;
        m_mesh.getTriangles(tcounts, tvers);
        
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
        B.begin_surface(pts.length(), tvers.length()/3);
        // vertices
        //typedef typename HDS::Vertex::Point Vertex;
        //typedef typename Vertex Point;
        typedef typename HDS::Vertex Vertex;
        typedef typename Vertex::Point Point;
        for ( unsigned i = 0 ; i < pts.length() ; i++ ) {
            
            /* original w/o IDs
            //HDS::Vertex_handle vh = B.add_vertex(Point(pts[i].x,pts[i].y,pts[i].z));
            //typename HDS::Vertex_handle vh = B.add_vertex(Point(pts[i].x,pts[i].y,pts[i].z));
            //vh->id = i;
            * end original */
            
            B.add_vertex(Point(pts[i].x,pts[i].y,pts[i].z));
        }
        
        // triangles
        for ( unsigned i = 0 ; i < tvers.length()/3 ; i++ ) {
            //HDS::Face_handle fh = B.begin_facet();
            //typename HDS::Face_handle fh = B.begin_facet();                   //original w/o IDs
            B.begin_facet();
            B.add_vertex_to_facet(tvers[3*i]);
            B.add_vertex_to_facet(tvers[3*i+1]);
            B.add_vertex_to_facet(tvers[3*i+2]);
            B.end_facet();
            //fh->id = i;
        }
        
        B.end_surface();
    }
    
private:
    MFnMesh m_mesh;
};

//why static?
static void mesh2polyhedron(MFnMesh &mesh, Polyhedron &P)
{
    Mesh_to_polyhedron<Polyhedron::HalfedgeDS> builder(mesh);
    P.delegate(builder);
    
    ///* Test for correctly built item
    CGAL::set_pretty_mode(cout);
    
    int vert_num = 0;
    cout << "created the polyhedron cgalPolyH:" << endl;
    for(Polyhedron::Vertex_iterator iter = P.vertices_begin(); iter != P.vertices_end(); ++iter) {
        vert_num++;
        //std::cout << iter->point() << std::endl;
        //cout << endl;
    }
    std::cout << "Number of CGAL vertices: " << vert_num << std::endl;
    
    int edge_num = 0;
    for(Polyhedron::Edge_iterator iter = P.edges_begin(); iter != P.edges_end(); ++iter) {
        edge_num++;
        //std::cout << iter->Vertex << std::endl;
        //cout << endl;
    }
    std::cout << "Number of CGAL edges: " << edge_num << std::endl;
    
    int half_num = 0;
    for(Polyhedron::Halfedge_iterator iter = P.halfedges_begin(); iter != P.halfedges_end(); ++iter) {
        half_num++;
        //std::cout << iter->Vertex << std::endl;
        //cout << endl;
    }
    std::cout << "Number of CGAL half-edges: " << half_num << std::endl;
     
    int facet_num = 0;
    for(Polyhedron::Facet_iterator iter = P.facets_begin(); iter != P.facets_end(); ++iter) {
        facet_num++;
        //std::cout << iter->Vertex << std::endl;
        //cout << endl;
    }
    std::cout << "Number of CGAL facets: " << facet_num << std::endl;
     //*/
}

//helper function for computing cotangent
static double cotan(double i) { return (1 / tan(i)); }

#endif

