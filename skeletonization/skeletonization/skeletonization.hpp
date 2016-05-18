//
//  skeletonization.hpp
//  skeletonization
//
//  Created by Gary Tse on 5/4/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef skeletonization_hpp
#define skeletonization_hpp

#include <stdio.h>
#include <iostream>

#include <maya/MPxCommand.h>
#include <maya/MGlobal.h>
#include <maya/MString.h>
#include <maya/MDagPath.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnTransform.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MFnMesh.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItMeshVertex.h>
#include <maya/MFloatArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MIOStream.h>
#include <maya/MFnPlugin.h>
#include <maya/MDagModifier.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Surface_mesh_deformation.h>

#include <CGAL/Eigen_solver_traits.h>
#include <Eigen/SparseCholesky>
#include <boost/graph/graph_traits.hpp>

#include "convertCGAL.cpp"


class Skeletonize:public MPxCommand {
public:
    
    //Maya functions
    Skeletonize();                                      // Constructor
    virtual ~Skeletonize();                             // Destructor
    virtual MStatus doIt(const MArgList &argList);
    //virtual MStatus redoIt();
    //virtual MStatus undoIt();
    //bool isUndoable() const;
    static void* creator();
    
    //CGAL functions
    typedef CGAL::Eigen_solver_traits<>     Eigen_Solver;
    typedef Eigen_Solver::Matrix            Eigen_matrix;
    typedef Eigen_Solver::Vector            Eigen_vector;

    void createLaplacian(Eigen_matrix &L, Polyhedron P);
    double calculate_volume(Polyhedron &cgalPolyH);
    
    void contract_geometry(Polyhedron &cgalPolyH, Eigen_matrix &L);
    void connectivity_surgery();
    void curve_refinement();
    
    //use for linear least square solver - part of contract_geometry
    //typedef typename CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> SolverTraits;
    
    //typedef typename CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<Eigen::MatrixXd>> SolverMatrix;
    //typedef CGAL::Eigen_matrix<Eigen::MatrixXd> MatrixXd;
    
private:
    MDagPath meshDagPath;
    MDagPath curveDagPath;
    
    Polyhedron cgalPolyH;
    
    /* May be able to use 'Matrix.set_coef(...)' instead
    typedef Eigen::Triplet<double> Triplet;
    mutable std::vector<Triplet> m_triplets;
    */
    
    /*
    MObject  	createWithEditPoints( const MPointArray &editPoints,
                                     unsigned int degree,
                                     Form agForm,
                                     bool create2D,
                                     bool createRational,
                                     bool uniformParam,
                                     MObject & parentOrOwner = MObject::kNullObj,
                                     MStatus* ReturnStatus = NULL );
     */
    MPointArray cEditPts;
};

#endif /* skeletonization_hpp */