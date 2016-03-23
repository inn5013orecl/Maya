//
//  doHelix.h
//  curveTest
//
//  Created by Gary Tse on 2/22/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef doHelix_h
#define doHelix_h

#include <math.h>
#include <maya/MSimple.h>
#include <maya/MIOStream.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MPointArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MFnPlugin.h>

class DoHelix:public MPxCommand {

public:
    DoHelix();                                          // Constructor
    ~DoHelix();                                         // Destructor
    virtual MStatus doIt(const MArgList& argList);      // function call
    static void* creator();                             // create instance of object
};

#endif /* doHelix_h */