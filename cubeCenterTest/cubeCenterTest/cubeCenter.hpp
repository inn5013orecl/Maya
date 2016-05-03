//
//  cubeCenter.hpp
//  cubeCenterTest
//
//  Created by Gary Tse on 4/19/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef cubeCenter_hpp
#define cubeCenter_hpp

#include <stdio.h>
#include <typeinfo>
#include <maya/MSimple.h>
#include <maya/MGlobal.h>
#include <maya/MString.h>
#include <maya/MDagPath.h>
#include <maya/MFnDagNode.h>
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

class CubeCenter:public MPxCommand {
public:
    CubeCenter();
    virtual ~CubeCenter();
    virtual MStatus doIt(const MArgList &argList);
    static void* creator();
};

#endif /* cubeCenter_hpp */
