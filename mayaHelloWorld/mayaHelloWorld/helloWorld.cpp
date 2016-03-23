//
//  helloWorld.cpp
//  mayaHelloWorld
//
//  Created by Gary Tse on 2/22/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

//#include <maya/MFnPlugin.h>
#include "helloWorld.hpp"

void* HelloWorld::creator() { return new HelloWorld; }

MStatus HelloWorld::doIt(const MArgList& argList) {
    MGlobal::displayInfo("Hello World!");
    return MS::kSuccess;
}