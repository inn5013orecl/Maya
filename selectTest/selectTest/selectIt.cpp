//
//  select.cpp
//  selectTest
//
//  Created by Gary Tse on 3/22/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef selectIt_cpp
#define selectIt_cpp

#include "selectIt.hpp"

void* SelectIt::creator() {return new SelectIt;}

MStatus SelectIt::doIt(const MArgList &argList) {
    
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
        
        //at the moment, only gets # of triangles for selected faces
        MItMeshPolygon triIter(meshDagPath);
        triIter.numTriangles(numTri);
    }
    
    /*
    for(MItMeshVertex vertexIter(thePolygon, &stat); !vertexIter.isDone(); vertexIter.next() ){
        verticesList.append(vertexIter.position());
    }*/
    
    cout << "Vertex List length: " << verticesList.length() << endl;
    
    for (int i = 0; i < verticesList.length(); i++) {
        cout << verticesList[i] << endl;
    }
    
    cout << "\nNumber of Triangles: " << numTri << endl;
    
    if (MS::kSuccess != stat) {
        cout << "Error in code." << endl;
    }
    else {
        cout << "Done." << endl;
    }
    return MS::kSuccess;
}

SelectIt::SelectIt() {}

SelectIt::~SelectIt() {}

MStatus
initializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin, "Gary Tse", "1.0.0", "Any");
    
    // Register elements here
    status = pluginFn.registerCommand("selectIter", SelectIt::creator);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return status;
}

MStatus
uninitializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin);
    
    // Deregister elements here
    status = pluginFn.deregisterCommand("selectIter");
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return status;
}

#endif