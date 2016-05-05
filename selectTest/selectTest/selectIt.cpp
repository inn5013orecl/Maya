//
//  select.cpp
//  selectTest
//
//  Created by Gary Tse on 3/22/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

// note on libc++ and Maya: cout does not work on Maya objects

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
    
    //cout << "\nNumber of Triangles: " << numTri << endl;
    
    /** Diagonals and node positioning **/
    
    // get diagonals to calculate center
    /*
    MPoint  first = verticesList[0]; //bottom left
    MPoint  second = verticesList[5]; //top right
    MPoint  diag = MPoint( (first.x + second.x)/2 , (first.y + second.y)/2 , (first.z + second.z)/2, 1 );
    
    //cout << "middle of cube: " << diag << endl;
    
    MObject centerPoint = MObject::kNullObj;
    MFnTransform centerPt;
    
    centerPoint = centerPt.create(MObject::kNullObj, &stat);
    //cout << "centerpt: " << centerPt.getTranslation(MSpace::kWorld) << endl;
    
    //cout << "diag vector: " << MVector(diag) << endl;
    
    double a = 0 , b= 2, c = 1;
    
    MVector diagonal = MVector(a,b,c);
    
    stat = centerPt.setTranslation(diagonal, MSpace::kWorld);
    
    if (MS::kSuccess != stat) {
        cout << "Error in setTrans." << endl;
    }
    
    //cout << "new centerpt: " << centerPt.getTranslation(MSpace::kWorld) << endl;
    
    
    MDagModifier cmd;
    
    MObject objLoc = cmd.createNode(MString("locator"),
                                    MObject::kNullObj,
                                    &stat);
    
    stat = cmd.doIt();
    
    //cout << objLoc.apiType() << endl; //110 = kTransform
    */
    /** End diagonals and node positioning **/
    
    if (MS::kSuccess != stat) {
        //cout << "Error in code." << endl;
    }
    else {
        //cout << "Done." << endl;
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