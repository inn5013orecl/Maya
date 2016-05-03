//
//  yTwist.hpp
//  yTwistDeform
//
//  Created by Gary Tse on 4/12/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef yTwist_h
#define yTwist_h

#include <stdio.h>
#include <string.h>
#include <maya/MIOStream.h>
#include <math.h>

#include <maya/MPxGeometryFilter.h>
#include <maya/MItGeometry.h>

#include <maya/MTypeId.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnPlugin.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MPoint.h>
#include <maya/MMatrix.h>

#define McheckErr(stat,msg)     \
    if ( MS::kSuccess != stat ) {   \
        cerr << msg;                \
        return MS::kFailure;        \
    }

class yTwist:public MPxGeometryFilter {
public:
    yTwist();
    virtual             ~yTwist();
    static  void*       creator();
    static  MStatus     initialize();
    // deformation function
    virtual MStatus deform(MDataBlock&      block,
                           MItGeometry&     iter,
                           const MMatrix&   mat,
                           unsigned int     multiIndex);

    
    // yTwist attributes
    static MObject angle;   // angle to twist
    
    static MTypeId id;
private:
};

#endif /* yTwist_h */
