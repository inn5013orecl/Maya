//
//  basicDepNode.hpp
//  basicDepNode
//
//  Created by Gary Tse on 3/7/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef basicDepNode_h
#define basicDepNode_h

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <maya/MGlobal.h>
#include <maya/MString.h>
#include <maya/MPxNode.h>
#include <maya/MTypeId.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnNumericAttribute.h>

class Sine: MPxNode {
public:
                        Sine();
    virtual             ~Sine();
    virtual MStatus     compute(const MPlug& plug, MDataBlock& data);
    static  void*       creator();
    static  MStatus     initialize();
    
    static  MObject     input;
    static  MObject     output;
    static  MTypeId     id;
};

#endif /* basicDepNode_h */
