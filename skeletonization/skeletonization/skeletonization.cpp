//
//  skeletonization.cpp
//  skeletonization
//
//  Created by Gary Tse on 5/4/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#include "convertCGAL.cpp"
#include "skeletonization.hpp"

void* Skeletonize::creator() { return new Skeletonize; }

Skeletonize::Skeletonize() {}
Skeletonize::~Skeletonize() {}

MStatus Skeletonize::doIt(const MArgList &argList) {
    
    Polyhedron cgalPolyH;
    
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
        
        //::mesh2polyhedron(meshFn, cgalPolyH);
        mesh2polyhedron(meshFn, cgalPolyH);
        
        //at the moment, only gets # of triangles for selected faces
        MItMeshPolygon triIter(meshDagPath);
        triIter.numTriangles(numTri);
    }
    
    /*
     for(MItMeshVertex vertexIter(thePolygon, &stat); !vertexIter.isDone(); vertexIter.next() ){
     verticesList.append(vertexIter.position());
     }*/
    
    cout << "Vertex List length: " << verticesList.length() << endl;
    
    std::cout << "Test string" << endl;
    
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

    
    if (MS::kSuccess != stat) {
        //cout << "Error in code." << endl;
    }
    else {
        //cout << "Done." << endl;
    }
    return MS::kSuccess;
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
