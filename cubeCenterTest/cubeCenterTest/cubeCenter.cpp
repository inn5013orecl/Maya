
//
//  cubeCenter.cpp
//  cubeCenterTest
//
//  Created by Gary Tse on 4/19/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#include "cubeCenter.hpp"

void* CubeCenter::creator() { return new CubeCenter; }

MStatus CubeCenter::doIt(const MArgList &argList) {

    MStatus stat;
    MDagPath meshDagPath;
    MFnDagNode nodeFn;
    MPointArray verticesList;
    MSelectionList list;
    MFnMesh meshFn;
    
    
    
    return MS::kSuccess;
}

MStatus
initializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin, "Gary Tse", "1.0.0", "Any");
    
    // Register elements here
    
    return status;
}

MStatus
uninitializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin);
    
    // Deregister elements here
    
    return status;
}