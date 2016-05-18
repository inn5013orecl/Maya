//
//  skeletonization.cpp
//  skeletonization
//
//  Created by Gary Tse on 5/4/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

//#include "convertCGAL.cpp"
#include "skeletonization.hpp"

void* Skeletonize::creator() { return new Skeletonize; }

Skeletonize::Skeletonize() {}
Skeletonize::~Skeletonize() {}

MStatus Skeletonize::doIt(const MArgList &argList) {
    
    MStatus stat;
    MDagPath meshDagPath;
    MFnDagNode nodeFn;
    //MObject thePolygon;
    MPointArray verticesList;
    MSelectionList list;
    MFnMesh meshFn;
    MIntArray triangles, trianglesIDs;
    int numTri;
    
    MGlobal::getActiveSelectionList(list);
    
    for (MItSelectionList listIter(list, MFn::kMesh, &stat); !listIter.isDone(); listIter.next()) {
        listIter.getDagPath(meshDagPath);
        meshFn.setObject(meshDagPath);
        
        meshFn.getPoints(verticesList, MSpace::kWorld);
        
        //convert Maya mesh to CGAL Polyhedron_3
        mesh2polyhedron(meshFn, cgalPolyH);
        
        //at the moment, only gets # of triangles for selected faces
        MItMeshPolygon triIter(meshDagPath);
        triIter.numTriangles(numTri);
    }
    
    // add ID's the Polyhedron
    CGAL::set_halfedgeds_items_id(cgalPolyH);
    
    // save original Polyhedron just in case
    Polyhedron cgalPolyH_original = cgalPolyH;
    
    std::size_t nvert = 0 /*, nrow = 0, ncol = 0 */;
    nvert = cgalPolyH.size_of_vertices();
    
    // initialize Laplacian Matrix (nvert x nvert)
    // all elements initialized to 0
    Eigen_matrix L(nvert,nvert);
    
    //contract geometry
    contract_geometry(cgalPolyH, L);
    
    /*
     for(MItMeshVertex vertexIter(thePolygon, &stat); !vertexIter.isDone(); vertexIter.next() ){
     verticesList.append(vertexIter.position());
     }*/
    
    cout << "Maya Vertex count: " << verticesList.length() << endl;
    
    /*
    for (int i = 0; i < verticesList.length(); i++) {
        //cout << verticesList[i] << endl;
        double vertList[4];
        verticesList[i].get(vertList);
        for (int j = 0; j < 4; j++) {
            cout << vertList[j] << ", ";
        }
        cout << endl;
        //cout << " -- "<< typeid(vertList[0]).name() << endl;
        //cout << endl;
    }
     */

    
    if (MS::kSuccess != stat) {
        //cout << "Error in code." << endl;
    }
    else {
        //cout << "Done." << endl;
    }
    return MS::kSuccess;
}

void Skeletonize::createLaplacian(Eigen_matrix &L, Polyhedron P) {
    
    std::size_t i, j;
    
    /***** i,j is edge in Polyhedron *****/
    // iterate through all halfedges
    for (Polyhedron::Halfedge_iterator iter = P.halfedges_begin(); iter != P.halfedges_end(); ++iter) {
        // grab vertex id ( i = end of halfedge)
        i = iter->opposite()->vertex()->id();               // source vertex ID
        j = iter->vertex()->id();                           // end vertex ID
        
        // grab facet incident to halfedge to calculate opposite angles
        Polyhedron::Halfedge_around_facet_circulator hc = iter->facet()->facet_begin();
        Polyhedron::Point_3 p, q, r;
        
        //make sure it's a triangle
        CGAL_assertion(hc->is_triangle());
        int stop = 0;
        do {
            if (hc->vertex() == iter->vertex()) {
                p = hc->vertex()->point();
                r = hc->opposite()->vertex()->point();
                q = hc->next()->vertex()->point();          // since circulator is counter-clockwise
                stop = 1;                                   // just in case break doesn't work for some reason
                break;
            }
        }while(++hc != iter->facet()->facet_begin() || stop == 1);
        
        // get angle
        CGAL::Cartesian<double>::Vector_3 v1(q, p),         // v1 = p - q
                                          v2(q, r);         // v2 = r - q
        
        // make unit vectors
        v1 = v1 / std::sqrt(v1 * v1);
        v2 = v2 / std::sqrt(v2 * v2);
        
        double angle1 = std::acos(v1 * v2);
        
        // do same for opposite facet
        Polyhedron::Halfedge_around_facet_circulator hc_o = iter->opposite()->facet()->facet_begin();
        
        //make sure it's a triangle
        CGAL_assertion(hc->is_triangle());
        int stop_o = 0;
        do {
            if (hc_o->vertex() == iter->opposite()->vertex()) {
                p = hc_o->vertex()->point();
                r = hc_o->opposite()->vertex()->point();
                q = hc_o->next()->vertex()->point();          // since circulator is counter-clockwise
                stop = 1;                                   // just in case break doesn't work for some reason
                break;
            }
        }while(++hc_o != iter->opposite()->facet()->facet_begin() || stop_o == 1);
        
        // get angle for opposite side
        CGAL::Cartesian<double>::Vector_3 v3(q, p),         // v3 = p - q
                                          v4(q, r);         // v4 = r - q
        
        // make unit vectors
        v3 = v3 / std::sqrt(v3 * v3);
        v4 = v4 / std::sqrt(v4 * v4);
        
        double angle2 = std::acos(v3 * v4);
        
        
        //do equation
        double w_ij = cotan(angle1) + cotan(angle2);
        
        // L.set(i,j, result_of_equation)
        L.set_coef(i, j, w_ij);
    }
    
    /***** diagonal of L matrix (L[i][i]) *****/
    // iterate through all vertices
    for (Polyhedron::Vertex_iterator iter = P.vertices_begin(); iter != P.vertices_end(); ++iter) {
        //term to be summed
        double w_ij = 0.0;
        
        // grab vertex id (i)
        i = iter->id();
        
        // iterate through all incident halfedges (note: vertex will be end of halfedge)
        Polyhedron::Halfedge_around_vertex_circulator hc = iter->vertex_begin();
        Polyhedron::Point_3 p, q, r, s, t, v;
        
        do {
            p = hc->vertex()->point();
            r = hc->opposite()->vertex()->point();
            q = hc->next()->vertex()->point();
            
            CGAL::Cartesian<double>::Vector_3 v1(q, p),         // v1 = p - q
                                              v2(q, r);         // v2 = r - q
            
            // make unit vectors
            v1 = v1 / std::sqrt(v1 * v1);
            v2 = v2 / std::sqrt(v2 * v2);
            
            double angle1 = std::acos(v1 * v2);
            
            s = hc->opposite()->vertex()->point();
            t = hc->opposite()->next()->vertex()->point();
            v = hc->vertex()->point();
            
            // get angle for opposite side
            CGAL::Cartesian<double>::Vector_3 v3(q, p),         // v3 = p - q
                                              v4(q, r);         // v4 = r - q
            
            // make unit vectors
            v3 = v3 / std::sqrt(v3 * v3);
            v4 = v4 / std::sqrt(v4 * v4);
            
            double angle2 = std::acos(v3 * v4);
            
            //do equation
            w_ij -= cotan(angle1) + cotan(angle2);
            
        }while(++hc != iter->vertex_begin());

        // L.set(i,i, result_of_equation)
        L.set_coef(i, i, w_ij);
        
    }
    
}

double Skeletonize::calculate_volume(Polyhedron &cgalPolyH) {
    return 0.0;
}

void Skeletonize::contract_geometry(Polyhedron &cgalPolyH, Eigen_matrix &L) {
    
    // initial values
    int iter_count = 0;
    double original_volume, current_volume, volume_ratio;
    double threshold = 1e-6;
    original_volume = calculate_volume(cgalPolyH);
    
    // create Laplacian matrix with new vertex positions every iteration
    do {
        createLaplacian(L, cgalPolyH);
        current_volume = calculate_volume(cgalPolyH);
        volume_ratio = current_volume / original_volume;
        ++iter_count;
    }while(volume_ratio >= threshold);
    
    std::cout << "\n===========\nContraction Complete\n===========\n\n";
}

void Skeletonize::connectivity_surgery() {
    std::cout << "\n===========\nConnectivity Surgery Complete\n===========\n\n";
}

void Skeletonize::curve_refinement() {
    std::cout << "\n===========\nCurve Refinement Complete\n===========\n\n";
}







MStatus
initializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin, "Gary Tse", "1.0.0", "Any");
    
    // Register elements here
    status = pluginFn.registerCommand("skeletonize", Skeletonize::creator);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return status;
}

MStatus
uninitializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin);
    
    // Deregister elements here
    status = pluginFn.deregisterCommand("skeletonize");
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return status;
}
