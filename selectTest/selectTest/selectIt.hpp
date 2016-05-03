//
//  select.hpp
//  selectTest
//
//  Created by Gary Tse on 3/22/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef selectIt_hpp
#define selectIt_hpp

#include <cstdio>
#include <iostream>
#include <typeinfo>
#include <maya/MSimple.h>
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

#include <CGAL/Surface_mesh_deformation.h>

class SelectIt:public MPxCommand {
public:
    SelectIt();
    virtual ~SelectIt();
    virtual MStatus doIt(const MArgList &argList);
    static void* creator();
};

#endif /* select_hpp */
