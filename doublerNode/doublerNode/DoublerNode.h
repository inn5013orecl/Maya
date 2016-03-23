//
//  DoublerNode.hpp
//  doublerNode
//
//  Created by Gary Tse on 3/7/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef DoublerNode_h
#define DoublerNode_h

#include <stdio.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MStatus.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MPxNode.h>

class DoublerNode : public MPxNode {
public:
    DoublerNode();
    virtual ~DoublerNode();
    virtual MStatus compute(const MPlug& plug, MDataBlock& data);
    static void* creator();
    static MStatus initialize();
    
    static MTypeId id;
    static MObject output;
    static MObject input;
    
};

#endif /* DoublerNode_h */
