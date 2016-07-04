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
    Polyhedron cgalPolyH_original = cgalPolyH;
    
    // initialize Laplacian Matrix (nvert x nvert)
    // all elements initialized to 0
    Eigen_matrix L(nvert,nvert);
    for (std::size_t i = 0; i < nvert; ++i) {
        for (std::size_t j = 0; j < nvert; ++j) {
            L.set(i, j, 0.0);                               //force everything to initilize to 0 due to random numbering errors in initialization
        }
    }
    
    //std::cout << "Laplacian Matrix L:\n" << L << std::endl;
    std::cout << "++++Finished initializing Laplacian matrix++++" << std::endl;
    
    //contract geometry
    std::cout << "++++BEFORE contract_geometry++++" << std::endl;
    contract_geometry(cgalPolyH, L, nvert, nface);
    
    /*
     for(MItMeshVertex vertexIter(thePolygon, &stat); !vertexIter.isDone(); vertexIter.next() ){
     verticesList.append(vertexIter.position());
     }*/
    
    cout << "Maya Vertex count: " << verticesList.length() << endl;
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
        L.set(i, j, w_ij);
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
        L.set(i, i, w_ij);
        
    }
    std::cout << "Laplacian Matrix L:\n" << L << std::endl;
    std::cout << "++++Finished Calculating Laplacian Matrix++++" << std::endl;
    
}

Polyhedron Skeletonize::quadratic_solver(Eigen_matrix W_l, Eigen_matrix W_h, Eigen_matrix L, Polyhedron P, std::size_t nvert) {

    Program qp(CGAL::EQUAL, false, 0, false, 0);
    
    // define A for A * x = b (quadratic programming constraint)
    // W_l * L
    Eigen::MatrixXd a_top = W_l * L;
    Eigen::MatrixXd a_bot = W_h;
    
    /*
    std::cout << "3rd column of first row of a_top: " << a_top.row(0)[2] << std::endl;
    std::cout << "4th column of 3rd row of a_top: " << a_top.row(2)[3] << std::endl;
    std::cout << "a_top:\n" << a_top << std::endl;
    std::cout << "before creating a_top vectors\n";
    */
    
    //NOTE: Use (the same) 1 row of a_top for every 3 rows of A, use next row of a_top for next 3 rows of A, repeat till no more rows of a_top
    
    Eigen::MatrixXd A(nvert * 6, nvert * 3);
    A.setZero();
    //std::cout << "W_L * L:\n" << a_top << std::endl;
    
    // top half of A
    std::size_t k = 0;
    for (std::size_t i = 0; i < 3 * nvert; ++i) {                           // row
        std::size_t m = 0;
        for (std::size_t j = 0; j < 3 * nvert; ++j) {                       // column
            if (j % 3 == i % 3) {
                qp.set_a(j, i, a_top.row(k)[m]);
                A(i,j) = a_top.row(k)[m];
                m++;
            }
            else {
                qp.set_a(j, i, 0);
            }
        }
        // next a_top row every 3 rows of A
        if ((i+1) % 3 == 0) k++;
    }
    
    // bottom half of A
    k = 0;
    for (std::size_t i = 3 * nvert; i < 6 * nvert; ++i) {                   // row
        std::size_t m = 0;
        for (std::size_t j = 0; j < 3 * nvert; ++j) {                       // column
            if (j % 3 == i % 3) {
                qp.set_a(j, i, a_bot.row(k)[m]);
                A(i,j) = a_bot.row(k)[m];
                m++;
            }
            else {
                qp.set_a(j, i, 0);
            }
        }
        // next a_top row every 3 rows of A
        if ((i+1) % 3 == 0) k++;
    }
    
    std::cout << "A Matrix:\n" << A << std::endl;
    
    // define b for A * x = b (quadratic programming constraint)
    Eigen::MatrixX3d vertices(nvert, 3);
    vertices.setZero();

    // create current iteration vertices matrix
    for (Polyhedron::Vertex_iterator iter = P.vertices_begin(); iter != P.vertices_end(); ++iter) {
        vertices(iter->id(), 0) = iter->point().x();
        vertices(iter->id(), 1) = iter->point().y();
        vertices(iter->id(), 2) = iter->point().z();
    }
    std::cout << "vertices matrix:\n" << vertices << std::endl;
    
    Eigen::VectorXd B(nvert * 6);
    Eigen::MatrixXd b_bot = W_h * vertices; //change L to current vertices matrix
    std::vector<double> b_bot_v;
    b_bot_v.clear();
    
    std::cout << "b_bot matrix:\n" << b_bot << std::endl;
    
    for (int i = 0, nRows = b_bot.rows(), nCols = b_bot.cols(); i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            b_bot_v.push_back(b_bot(i,j));
        }
    }
    
    for (std::size_t i = 0; i < 3 * nvert; ++i) {
        qp.set_b(i, 0);
        B(i) = 0;
    }
    k = 0;
    for (std::size_t i = 3 * nvert; i < 6 * nvert; ++i) {
        qp.set_b(i, b_bot_v[k]);
        B(i) = b_bot_v[k];
        k++;
    }
    
    std::cout << "B vector:\n" << B << std::endl;
    // solve for minimized energy
    //Eigen::MatrixXd W_h_squared(nvert,nvert);
    //W_h_squared = W_h * W_h;
    
    Eigen::MatrixXd values(nvert,nvert);
    //std::vector<double> c_values;

    //c_values.clear();
    values.setZero();
    
    // set D
    for (std::size_t i = 0; i < nvert; ++i) {                               // row
        for (std::size_t j = 0; j < nvert; ++j) {                           // column
            if (i == j) {
                values(i,j) += (W_h(i,j) * W_h(i,j)) + (a_top(i,j) * a_top(i,j));
            }
            else {
                values(i,j) += (a_top(i,j) * a_top(i,j));
            }
            //qp.set_d(j, i, 2 * value);
        }
    }
    
    for (std::size_t i = 0; i < nvert; ++i) {                               // row
        for (std::size_t j = 0; j < nvert; ++j) {                           // column
            qp.set_d(j, i, 2 * values(i,j));
        }
        
    }
    
    // set C
    /*
    for (std::size_t i = 0; i < nvert; ++i) {
        c_values.push_back(vertices(i,0));
        c_values.push_back(vertices(i,1));
        c_values.push_back(vertices(i,2));
    }
    
    for (std::size_t i = 0; i < nvert * 3; ++i) {
        qp.set_c(i, c_values[i]);
    }
    */
    
    for (std::size_t i = 0; i < nvert; i+=3) {
        qp.set_c(i,   vertices(i,0));
        qp.set_c(i+1, vertices(i,1));
        qp.set_c(i+2, vertices(i,2));
    }
    
    qp.set_c0(1);                                                   // c_0 needs to be sum of all of them; w_h^2 * (vertices(i).x^2 + vertices(i).y^2 + vertices(i).z^2)
    Polyhedron poly = P;
    
    Solution s = CGAL::solve_quadratic_program(qp, ET());
    auto solution_iter = s.variable_values_begin();
    s.variable_numerators_begin();
    //typedef CGAL::Cartesian<double>                 Kernel;
    //typedef Kernel::Point_3                         Point;
    for (Polyhedron::Vertex_iterator iter = poly.vertices_begin(); iter != poly.vertices_end(); ++iter) {
        if (solution_iter == s.variable_values_end()) {
            std::cout << "~~~~ERROR: solution values reached before polyhedron written~~~~" << std::endl;
            break;
        }
        
        double x = CGAL::to_double(solution_iter->numerator());
        double y = CGAL::to_double((solution_iter+1)->numerator())/CGAL::to_double((solution_iter+1)->denominator());
        double z = CGAL::to_double((solution_iter+2)->numerator())/CGAL::to_double((solution_iter+2)->denominator());
        iter->point() = Polyhedron::Point_3(x,y,z);
        solution_iter += 3;
    }

    return poly;
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
        //std::cout << "Initialized Matrix m:\n" << m << std::endl;
        
        // 1st row
        m.set(0, 0, p1.x());
        m.set(0, 1, p2.x());
        m.set(0, 2, p3.x());
        m.set(0, 3, 0);
        // 2nd row
        m.set(1, 0, p1.y());
        m.set(1, 1, p2.y());
        m.set(2, 2, p3.y());
        m.set(3, 3, 0);
        // 3rd row
        m.set(2, 0, p1.z());
        m.set(2, 1, p2.z());
        m.set(2, 2, p3.z());
        m.set(2, 3, 0);
        // 4th row
        m.set(3, 0, 1);
        m.set(3, 1, 1);
        m.set(3, 2, 1);
        m.set(3, 3, 1);
        
        //std::cout << "Post-set Matrix m:\n" << m << std::endl;
        
        double a = m.determinant();
        
        a = std::abs(a) / 6;
        
        // 6 * Volume = sum of determinants of (triangular) faces of polyhedron
        volume += a;

        // probably doesn't work due to volume calculations and triangle orientations
        // orientation is defined by location of CGAL::ORIGIN in relation to plane defined by p1,p2,p3
        //CGAL::Tetrahedron_3<CGAL::Cartesian<double>> tetra(p1,p2,p3,CGAL::ORIGIN);
        //tetra.volume();
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

void Skeletonize::contract_geometry(Polyhedron &cgalPolyH, Eigen_matrix &L, std::size_t nvert, std::size_t nface) {
    
    // initial values
    std::cout << "++++Begin Contract Geometry++++" << std::endl;
    int iter_count = 0;
    double original_volume = 0.0, current_volume = 0.0, volume_ratio = 0.0, current_or_area = 0.0;
    double threshold = 1e-6;
    std::cout << "++++BEFORE 1st calculate_volume++++" << std::endl;
    original_volume = calculate_volume(cgalPolyH);
    std::vector<double> original_or_area;
    original_or_area.clear();
    
    // calculate original one ring areas for every vertex;
    for (Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
        original_or_area.push_back(calculate_one_ring_area(iter));
    }
    
    std::cout << "Size (# of items in) of original_or_area: " << original_or_area.size() << std::endl;
    
    Eigen_matrix W_h(nvert, nvert), W_l(nvert, nvert), W_h_orig(nvert, nvert);

    // calculate initial avg face area term for W_l
    double f_area = 0.0;
    for (Polyhedron::Facet_iterator iter = cgalPolyH.facets_begin(); iter != cgalPolyH.facets_end(); ++iter) {
        CGAL_assertion(iter->is_triangle());
        Polyhedron::Point_3 p1,p2,p3;
        
        p1 = iter->halfedge()->vertex()->point();
        p2 = iter->halfedge()->next()->vertex()->point();
        p3 = iter->halfedge()->next()->next()->vertex()->point();
        
        CGAL::Triangle_3<CGAL::Cartesian<double>> t(p1,p2,p3);
        //std::cout << "current face's area: " << sqrt(t.squared_area()) << std::endl;
        f_area += sqrt(t.squared_area());
    }
    f_area = f_area / nface;
    
    for (std::size_t i = 0; i < nvert; ++i) {
        for (std::size_t j = 0; j < nvert; ++j) {
            if (i == j) {
                W_h.set(i, i, 1.0);
                W_l.set(i, i, pow(10, -3) * sqrt(f_area));
                // save original W_h for update purposes
                W_h_orig.set(i, i, 1.0);
            }
            else {
                W_h.set(i, j, 0.0);
                W_l.set(i, j, 0.0);
                // save original W_h for update purposes
                W_h_orig.set(i, j, 0.0);
            }
        }
    }
    
    std::cout << "W_h matrix:\n" << W_h << std::endl;
    std::cout << "W_l matrix:\n" << W_l << std::endl;
    
    
    //do {
        // create Laplacian matrix with new vertex positions every iteration
        createLaplacian(L, cgalPolyH);
    
        // Quadratic function minimization
        cgalPolyH = quadratic_solver(W_l, W_h, L, cgalPolyH, nvert);
    
        // update weights
        W_l *= 2.0;
        for (Polyhedron::Vertex_iterator iter = cgalPolyH.vertices_begin(); iter != cgalPolyH.vertices_end(); ++iter) {
            std::size_t i;
            i = iter->id();
            current_or_area = calculate_one_ring_area(iter);
            W_h.set(i,i, sqrt(original_or_area[i]/current_or_area));
            //W_h.set(i,i, W_h_orig(i,i) * sqrt(original_or_area[i]/current_or_area[i]));       // w_h_orig diagonals = 1.0
        }
    
        // calculate current volume after each iteration
        std::cout << "++++BEFORE 2nd calculate volume++++" << std::endl;
        current_volume = calculate_volume(cgalPolyH);
        std::cout << "++++AFTER 2nd calculate volume++++" << std::endl;
        volume_ratio = current_volume / original_volume;
        std::cout << "++++AFTER calculating volume RATIO++++" << std::endl;
        
        ++iter_count;
    //}while(volume_ratio >= threshold || iter_count >= 30);                                    // iteration count max 30 as a breakpoint
    
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
