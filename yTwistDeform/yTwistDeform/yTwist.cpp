//
//  yTwist.cpp
//  yTwistDeform
//
//  Created by Gary Tse on 4/12/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#include "yTwist.h"

MTypeId yTwist::id(0x00000002);

MObject yTwist::angle;

// constructor
yTwist::yTwist() {}

// destructor
yTwist::~yTwist() {}

// create plugin
void* yTwist::creator() { return new yTwist(); }

MStatus yTwist::initialize() {
    
    MFnNumericAttribute nAttr;
    angle = nAttr.create("angle", "ta", MFnNumericData::kDouble);
    
    nAttr.setDefault(0.0);
    nAttr.setKeyable(true);
    addAttribute(angle);
    
    attributeAffects(angle, outputGeom);
    
    return MS::kSuccess;
}

//twist points around Y axis
MStatus yTwist::deform(MDataBlock&      block,
                       MItGeometry&     iter,
                       const MMatrix&   /*mat*/,
                       unsigned int     /*multiIndex*/) {
    
    MStatus status = MS::kSuccess;
    
    MDataHandle angleData =  block.inputValue(angle, &status);
    McheckErr(status, "Error with angleData\n");
    double magnitude = angleData.asDouble();
    
    
    // envelopeData required for global scaling factors
    MDataHandle envolopeData = block.inputValue(envelope, &status);
    McheckErr(status, "Error with evelopeData\n");
    float env = envolopeData.asFloat();
    
    for(; !iter.isDone(); iter.next()) {
        MPoint pt = iter.position();
        
        double ang = magnitude * pt.y * env;
        if (ang != 0.0) {
            double c = cos(ang);
            double s = sin(ang);
            double t = pt.x * c - pt.z * s;

            //set point's X and Z axis
            pt.z = pt.x * s + pt.z * c;
            pt.x = t;
        }
        
        iter.setPosition(pt);
    }
    
    return status;
}

MStatus
initializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin, "Gary Tse", "3.0.0", "Any");
    // Register elements here
    status = pluginFn.registerNode("ytwist", yTwist::id, yTwist::creator, yTwist::initialize, MPxNode::kDeformerNode);
    
    return status;
}

MStatus
uninitializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin);
    // Deregister elements here
    status = pluginFn.deregisterNode(yTwist::id);
    
    return status;
}