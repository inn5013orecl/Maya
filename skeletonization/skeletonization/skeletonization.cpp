//
//  skeletonization_2.cpp
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
        //MItMeshPolygon triIter(meshDagPath);
        //triIter.numTriangles(numTri);
    }
    
    // add ID's the Polyhedron
    CGAL::set_halfedgeds_items_id(cgalPolyH);
    std::cout << "++++Finished adding IDs++++" << std::endl;
    
    std::size_t nvert = 0, nface = 0; /*, nrow = 0, ncol = 0 */;
    nvert = cgalPolyH.size_of_vertices();
    nface = cgalPolyH.size_of_facets();
    std::cout << "Number of vertices = " << nvert << std::endl;
    
    CGAL::set_ascii_mode(cout);
    
    /*
     for (Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
     std::cout << "vertex " << iter->id() << ": " << iter->point() << std::endl;
     }
     std::cout<<std::endl;
     
     std::size_t edge_count = 0;
     for (Polyhedron::Halfedge_iterator iter = cgalPolyH.halfedges_begin(); iter != cgalPolyH.halfedges_end(); ++iter) {
     std::cout << "Halfedge " << edge_count << ":\nsource: " << iter->opposite()->vertex()->point() << ", end: " << iter->vertex()->point() << std::endl;
     std::cout << "Halfedge next: " << iter->next()->vertex()->point() <<std::endl;
     edge_count++;
     }
     std::cout<<std::endl;
     */
    
    // save original Polyhedron just in case
    //Polyhedron cgalPolyH_original = cgalPolyH;
    
    // initialize Laplacian Matrix (nvert x nvert)
    // all elements initialized to 0
    Eigen_matrix L(nvert,nvert);
    L.setZero();
    
    //  //std::cout << "Laplacian Matrix L:\n" << L << std::endl;
    //std::cout << "++++Finished initializing Laplacian matrix++++" << std::endl;
    
    //contract geometry
    std::cout << "++++BEFORE contract_geometry++++" << std::endl;
    std::clock_t start = std::clock();
    
    contract_geometry(cgalPolyH, L, nvert, nface);
    
    double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    
    /*
     for(MItMeshVertex vertexIter(thePolygon, &stat); !vertexIter.isDone(); vertexIter.next() ){
     verticesList.append(vertexIter.position());
     }*/
    
    std::cout << "Maya Vertex count: " << verticesList.length() << std::endl;
    
    std::cout << "Time to contract: " << duration << " seconds" << std::endl;
    
    print_obj_file(cgalPolyH);
    //cout << "Maya vertices list:" << endl;
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
                q = hc->next()->vertex()->point();              // since circulator is counter-clockwise
                stop = 1;                                       // just in case break doesn't work for some reason
                break;
            }
        }while(++hc != iter->facet()->facet_begin() || stop == 1);
        
        // get angle
        CGAL::Cartesian<double>::Vector_3 v1(q, p),             // v1 = p - q
        v2(q, r);             // v2 = r - q
        
        // make unit vectors
        v1 = v1 / std::sqrt(v1.squared_length());
        v2 = v2 / std::sqrt(v2.squared_length());
        
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
                q = hc_o->next()->vertex()->point();            // since circulator is counter-clockwise
                stop = 1;                                       // just in case break doesn't work for some reason
                break;
            }
        }while(++hc_o != iter->opposite()->facet()->facet_begin() || stop_o == 1);
        
        // get angle for opposite side
        CGAL::Cartesian<double>::Vector_3 v3(q, p),             // v3 = p - q
        v4(q, r);             // v4 = r - q
        
        // make unit vectors
        v3 = v3 / std::sqrt(v3.squared_length());
        v4 = v4 / std::sqrt(v4.squared_length());
        
        double angle2 = std::acos(v3 * v4);
        
        
        //do equation
        double w_ij = cotan(angle1) + cotan(angle2);
        
        // L.set(i,j, result_of_equation)
        //L.set_coef(i, j, w_ij);
        L(i, j) = w_ij;
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
            
            
            //CGAL::Cartesian<double>::Vector_3
            CGAL::Vector_3<CGAL::Cartesian<double>> v1(q, p),   // v1 = p - q
            v2(q, r);   // v2 = r - q
            
            // make unit vectors
            v1 = v1 / std::sqrt(v1.squared_length());
            v2 = v2 / std::sqrt(v2.squared_length());
            
            double angle1 = std::acos(v1 * v2);
            
            s = hc->opposite()->vertex()->point();
            t = hc->opposite()->next()->vertex()->point();
            v = hc->vertex()->point();
            
            // get angle for opposite side
            CGAL::Vector_3<CGAL::Cartesian<double>> v3(q, p),   // v3 = p - q
            v4(q, r);   // v4 = r - q
            
            // make unit vectors
            v3 = v3 / std::sqrt(v3.squared_length());
            v4 = v4 / std::sqrt(v4.squared_length());
            
            double angle2 = std::acos(v3 * v4);
            
            //do equation
            w_ij -= cotan(angle1) + cotan(angle2);
            
        }while(++hc != iter->vertex_begin());
        
        // L.set(i,i, result_of_equation)
        //L.set_coef(i, i, w_ij);
        L(i, i) = w_ij;
        
    }
    //std::cout << "Laplacian Matrix L:\n" << L << std::endl;
    std::cout << "++++Finished Calculating Laplacian Matrix++++" << std::endl;
    
}

void Skeletonize::lls_solver(Eigen_matrix W_l, Eigen_matrix W_h, Eigen_matrix L, Polyhedron &cgalPolyH, std::size_t nvert, std::vector<bool> degenerate_vertex) {
    
    // create A
    Eigen_matrix a_top = W_l * L;
    //std::cout << "Before Matrix A" << std::endl;
    SolverTraits::Matrix A(nvert * 2, nvert);
    std::cout << "Matrix A created" << std::endl;
    
    for (Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
        Polyhedron::Halfedge_around_vertex_circulator hc = iter->vertex_begin();
        std::size_t i = iter->id();
        if (degenerate_vertex[iter->id()]) {
            // bottom half
            A.set_coef(i + nvert, i, 1/ 0.0000001);
            // top half
            do {
                A.set_coef(i, hc->opposite()->vertex()->id(), 0);
            }while (++hc != iter->vertex_begin());
        }
        else {
            // bottom half
            A.set_coef(i + nvert, i, W_h(i,i));
            // top half
            do {
                A.set_coef(i, hc->opposite()->vertex()->id(), a_top(i,hc->opposite()->vertex()->id()));
            }while (++hc != iter->vertex_begin());
        }
        A.set_coef(i, i, a_top(i,i));
    }

    /* old code
    for (std::size_t i = 0; i < nvert; ++i) {
        for (std::size_t j = 0; j < nvert; ++j) {
            // top half
            A.set_coef(i, j, a_top(i,j));
            // bottom half
            A.set_coef(i + nvert, j, W_h(i,j));
        }
    }
    */
    
    std::cout << "Finished setting A" << std::endl;
    
    // create B
    SolverTraits::Vector X(nvert), Bx(nvert * 2);
    SolverTraits::Vector Y(nvert), By(nvert * 2);
    SolverTraits::Vector Z(nvert), Bz(nvert * 2);
    
    // create current iteration vertices matrix
    Eigen_matrix vertices(nvert, 3);
    vertices.setZero();
    
    for (Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
        vertices(iter->id(), 0) = iter->point().x();
        vertices(iter->id(), 1) = iter->point().y();
        vertices(iter->id(), 2) = iter->point().z();
    }
    //std::cout << "Pre-solving vertices matrix:\n" << vertices << std::endl;
    
    for (std::size_t i = 0; i < nvert; ++i) {
        // top half
        Bx[i] = 0;
        By[i] = 0;
        Bz[i] = 0;
        
        // bottom half
        if (degenerate_vertex[i]) {
            Bx[i + nvert] = (1/0.0000001) * vertices(i,0);
            By[i + nvert] = (1/0.0000001) * vertices(i,1);
            Bz[i + nvert] = (1/0.0000001) * vertices(i,2);
        }
        else {
            Bx[i + nvert] = W_h(i,i) * vertices(i,0);
            By[i + nvert] = W_h(i,i) * vertices(i,1);
            Bz[i + nvert] = W_h(i,i) * vertices(i,2);
        }
    }
    
    // solve At * A * X = At * B
    m_solver.normal_equation_factor(A);
    m_solver.normal_equation_solver(Bx, X);
    m_solver.normal_equation_solver(By, Y);
    m_solver.normal_equation_solver(Bz, Z);
    
    // copy into original mesh
    std::vector<Polyhedron::Point_3> pts;
    pts.reserve(nvert);
    
    for (std::size_t i = 0; i < nvert; ++i) {
        pts.push_back(Polyhedron::Point_3(X[i],Y[i],Z[i]));
    }
    
    std::copy(pts.begin(), pts.end(), cgalPolyH.points_begin());
    std::cout << "+++++Linear Least Squared Solver completed+++++" << std::endl;
    
    CGAL_postcondition(cgalPolyH.is_valid());
    
}

double Skeletonize::calculate_volume(Polyhedron P) {
    
    double volume = 0.0;
    
    for (Polyhedron::Facet_iterator iter = P.facets_begin(); iter != P.facets_end(); ++iter) {
        //make sure it's a triangle
        CGAL_assertion(iter->is_triangle());
        Polyhedron::Point_3 p1, p2, p3;
        
        p1 = iter->halfedge()->vertex()->point();
        p2 = iter->halfedge()->next()->vertex()->point();
        p3 = iter->halfedge()->next()->next()->vertex()->point();
        
        //CGAL::Linear_algebraCd<double>::Matrix m(4);
        //double a = CGAL::Linear_algebraCd<double>::determinant(m);
        
        Eigen_matrix m(4, 4);
        //  //std::cout << "Initialized Matrix m:\n" << m << std::endl;
        
        // 1st row
        m(0, 0) = p1.x();
        m(0, 1) = p2.x();
        m(0, 2) = p3.x();
        m(0, 3) = 0;
        // 2nd row
        m(1, 0) = p1.y();
        m(1, 1) = p2.y();
        m(1, 2) = p3.y();
        m(1, 3) = 0;
        // 3rd row
        m(2, 0) = p1.z();
        m(2, 1) = p2.z();
        m(2, 2) = p3.z();
        m(2, 3) = 0;
        // 4th row
        m(3, 0) = 1;
        m(3, 1) = 1;
        m(3, 2) = 1;
        m(3, 3) = 1;
        
        //  //std::cout << "Post-set Matrix m:\n" << m << std::endl;
        
        double a = m.determinant();
        
        a = std::abs(a) / 6;
        //a = a / 6;
        
        // 6 * Volume = sum of determinants of (triangular) faces of polyhedron
        volume += a;
        
    }
    std::cout << "volume: " << volume << std::endl;
    return volume;
}

double Skeletonize::calculate_one_ring_area(Polyhedron::Vertex_iterator vi) {
    
    double or_area = 0.0;
    Polyhedron::Halfedge_around_vertex_const_circulator hc = vi->vertex_begin();
    
    do {
        Polyhedron::Point_3 a,b,c;
        a = hc->vertex()->point();
        b = hc->next()->vertex()->point();
        c = hc->opposite()->vertex()->point();
        
        CGAL::Triangle_3<CGAL::Cartesian<double>> t(a,b,c);
        or_area += sqrt(t.squared_area());
    }while(++hc != vi->vertex_begin());
    
    return or_area;
}

bool Skeletonize::detect_degeneracies(Polyhedron::Vertex_iterator vertex, double min_edge_length) {
    
    std::set<Polyhedron::Vertex_handle> vertices_in_disk;
    std::set<Polyhedron::Halfedge_handle> edges_in_disk;
    std::set<Polyhedron::Facet_handle> facets_in_disk;
    
    vertices_in_disk.clear(); //edges_in_disk.clear(); facets_in_disk.clear();

    std::queue<Polyhedron::Vertex_handle> Q;
    
    // put vertex in disk and set disk checks for halfedges & faces
    Q.push(vertex);
    vertices_in_disk.insert(vertex);
    //std::size_t ct = 0, dt = 0;
    while (!Q.empty()) {
        //std::cout << "Q not empty: " << ++ct << std::endl;
        Polyhedron::Vertex_handle v = Q.front();
        Q.pop();
        //dt = 0;
        Polyhedron::Halfedge_around_vertex_circulator hc = v->vertex_begin();
        do {
            //std::cout << "in do-while: " << dt++ << std::endl;
            Polyhedron::Vertex_handle new_v = hc->opposite()->vertex();
            if (!vertices_in_disk.count(new_v)) {
                double distance = std::sqrt(CGAL::squared_distance(vertex->point(), new_v->point()));
                if (distance < min_edge_length) {
                    Q.push(new_v);
                    vertices_in_disk.insert(new_v);
                }
            }
        } while (++hc != v->vertex_begin());
    }
    //std::cout << "Q is now empty for current vertex" << std::endl;
    
    // check disk elements for euler characteristics
    for (std::set<Polyhedron::Vertex_handle>::iterator vi = vertices_in_disk.begin(); vi != vertices_in_disk.end(); ++vi) {
        Polyhedron::Vertex_handle vert = *vi;
        
        Polyhedron::Halfedge_around_vertex_circulator hv = vert->vertex_begin();
        do {
            Polyhedron::Vertex_handle v = hv->opposite()->vertex();
            if (vertices_in_disk.find(v) != vertices_in_disk.end()) {
                edges_in_disk.insert(hv->opposite());
                edges_in_disk.insert(hv);
            }
            
            bool face_inside = true;
            Polyhedron::Halfedge_around_facet_circulator hf = hv->facet_begin();
            do {
                Polyhedron::Vertex_handle v = hf->vertex();
                if (vertices_in_disk.find(v) == vertices_in_disk.end()) {
                    face_inside = false;
                    break;
                }
            } while (++hf != hv->facet_begin());
            
            if (face_inside) {
                facets_in_disk.insert(hv->facet());
            }
            
        } while (++hv != vert->vertex_begin());
    }
    //std::cout << "Finished checking disk" << std::endl;
    
    std::size_t V = vertices_in_disk.size();
    std::size_t E = edges_in_disk.size() / 2;
    std::size_t F = facets_in_disk.size();
    std::size_t euler = V + F - E;
    if (euler != 1)
    {
        return true;
    }
    return false;
    
}

void Skeletonize::contract_geometry(Polyhedron &cgalPolyH, Eigen_matrix &L, std::size_t nvert, std::size_t nface) {
    
    // initial values
    std::cout << "++++Begin Contract Geometry++++" << std::endl;
    int iter_count = 0;
    double original_volume = 0.0, current_volume = 0.0, volume_ratio = 100000.0, current_or_area = 0.0;
    double threshold = 0.000001;
    double min_edge_length = init_min_edge_length();
    std::cout << "++++BEFORE 1st calculate_volume++++" << std::endl;
    original_volume = calculate_volume(cgalPolyH);
    std::vector<double> original_or_area;
    original_or_area.clear();
    std::vector<bool> degenerate_vertex(nvert,0);
    
    // calculate original one ring areas for every vertex;
    for (Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
        original_or_area.push_back(calculate_one_ring_area(iter));
    }
    
    std::cout << "Size (# of items in) of original_or_area: " << original_or_area.size() << std::endl;
    
    Eigen_matrix W_l(nvert, nvert), W_h(nvert, nvert), W_h_orig(nvert, nvert);
    
    // calculate initial avg face area term for W_l
    double f_area = 0.0;
    for (Polyhedron::Facet_iterator iter = cgalPolyH.facets_begin(); iter != cgalPolyH.facets_end(); ++iter) {
        CGAL_assertion(iter->is_triangle());
        Polyhedron::Point_3 p1,p2,p3;
        
        p1 = iter->halfedge()->vertex()->point();
        p2 = iter->halfedge()->next()->vertex()->point();
        p3 = iter->halfedge()->next()->next()->vertex()->point();
        
        CGAL::Triangle_3<CGAL::Cartesian<double>> t(p1,p2,p3);
        //  //std::cout << "current face's area: " << sqrt(t.squared_area()) << std::endl;
        f_area += sqrt(t.squared_area());
    }
    f_area = f_area / nface;
    
    W_l.setZero(); W_h.setZero(); W_h_orig.setZero();
    
    for (std::size_t i = 0; i < nvert; ++i) {
        for (std::size_t j = 0; j < nvert; ++j) {
            if (i == j) {
                W_h(i, i) = 1.0;
                W_l(i, i) = pow(10, -3) * sqrt(f_area);
                // save original W_h for update purposes
                W_h_orig(i, i) = 1.0;
            }
            // all non-diagonals are zero
        }
    }
    
    //  //std::cout << "W_h matrix:\n" << W_h << std::endl;
    //  //std::cout << "W_l matrix:\n" << W_l << std::endl;
    
    
    while (volume_ratio > threshold && iter_count < 25) {                                             // iteration count max 30 as a breakpoint
        // create Laplacian matrix with new vertex positions every iteration
        createLaplacian(L, cgalPolyH);
    
        // Solver for linear least squared
        lls_solver(W_l, W_h, L, cgalPolyH, nvert, degenerate_vertex);
    
        /*
         std::cout << "post solution vertices:" << std::endl;
        for(Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
            std::cout << iter->point() << std::endl;
            //cout << endl;
        }
        */
        
        // check for degeneracies
        for (Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
            degenerate_vertex[iter->id()] = detect_degeneracies(iter, min_edge_length);
        }
        std::cout << "+++++Done checking degeneracies++++" << std::endl;
        // update weights
        W_l *= 2.0;
        for (Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
            std::size_t i;
            i = iter->id();
            current_or_area = calculate_one_ring_area(iter);
            W_h(i,i) = sqrt(original_or_area[i]/current_or_area);                                       // w_h_orig diagonals = 1.0
            //W_h.set(i,i, W_h_orig(i,i) * sqrt(original_or_area[i]/current_or_area[i]));
        }
        
        // calculate current volume after each iteration
        //std::cout << "++++BEFORE 2nd calculate volume++++" << std::endl;
        current_volume = calculate_volume(cgalPolyH);
        //std::cout << "++++AFTER 2nd calculate volume++++" << std::endl;
        volume_ratio = current_volume / original_volume;
        std::cout << "volume ratio: " << volume_ratio << std::endl;
        //std::cout << "++++AFTER calculating volume RATIO++++" << std::endl;
        
        ++iter_count;
    }
    
    std::cout << "iter count: " << iter_count << std::endl;
    
    if (iter_count >= 30) {
        std::cout << "max iterations of 30 reached" << std::endl;
    }
    
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
