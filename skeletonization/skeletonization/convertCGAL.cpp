//
//  convertCGAL.cpp
//  skeletonization
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
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

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
typedef CGAL::Polyhedron_3<Kernel, My_items> Polyhedron;


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
            //HDS::Vertex_handle vh = B.add_vertex(Point(pts[i].x,pts[i].y,pts[i].z));
            typename HDS::Vertex_handle vh = B.add_vertex(Point(pts[i].x,pts[i].y,pts[i].z));
            vh->id = i;
        }
        
        // triangles
        for ( unsigned i = 0 ; i < tvers.length()/3 ; i++ ) {
            //HDS::Face_handle fh = B.begin_facet();
            typename HDS::Face_handle fh = B.begin_facet();
            B.add_vertex_to_facet(tvers[3*i]);
            B.add_vertex_to_facet(tvers[3*i+1]);
            B.add_vertex_to_facet(tvers[3*i+2]);
            B.end_facet();
            fh->id = i;
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
    
    CGAL::set_pretty_mode(cout);
    cout << "created the polyhedron cgalPolyH:" << endl;
    for(Polyhedron::Vertex_iterator iter = P.vertices_begin(); iter != P.vertices_end(); ++iter) {
        std::cout << iter->point() << std::endl;
        cout << endl;
    }
}

#endif

