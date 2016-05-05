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


class Skeletonize:public MPxCommand {
public:
    Skeletonize();
    virtual ~Skeletonize();
    virtual MStatus doIt(const MArgList &argList);
    //virtual MStatus redoIt();
    //virtual MStatus undoIt();
    //bool isUndoable() const;
    static void* creator();

private:
    MDagPath meshDagPath;
    MDagPath curveDagPath;
    
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